#include <stdio.h>
#include <stdlib.h>

#if defined(__linux__ ) || defined(__APPLE__)
    #include <pthread.h>
    #define THREAD_FUNC_DEF(func_name) void *func_name(void *arguments)
    #define CREATE_THREAD(func_name) pthread_create(&tid_array[i], NULL, \
                                                    func_name, \
                                                    &thread_args[i])
    #define RETURN_SECTION NULL
#elif defined(_WIN32)
    #include <windows.h>
    #define THREAD_FUNC_DEF(func_name) DWORD WINAPI func_name(LPVOID arguments)
    #define CREATE_THREAD(func_name) hThreadArray[i] = CreateThread(NULL, 0, \
                                                  func_name, &thread_args[i], \
                                                  0, &dwThreadIdArray[i])
    #define RETURN_SECTION 0
#endif

#include <limits.h>
#include <stdbool.h>

// Our includes
#include "buffers.c"
#include "parsers.c"

// How big the file chunking will be
#define RECORD_CHUNK 1024*8 // 1024*8 gives the best results on Guillem laptop
#define MAX_RING_BUF 4096

int c_fseek(FILE *filehandle, long int offset)
{
    return fseek(filehandle, offset, SEEK_SET);
}

static inline void load_buffer(uint32_t *pbuffer, FILE *fhandle)
{
    if(fread(pbuffer, RECORD_CHUNK, sizeof(uint32_t), fhandle)==0) {
        if (ferror(fhandle)){
            perror("Error detected while reading file.");
            exit(0);
        }
    }
}

static inline void check_and_grow_buf(ring_buf_t *cbuf, uint64_t timetag,
                                      uint64_t correlation_window) {
    if ( (timetag-ring_buf_oldest(cbuf)) < correlation_window && cbuf->count == cbuf->size && cbuf->size < MAX_RING_BUF) {
        ring_buf_grow(cbuf);
    }
}


static inline bool next_record(FILE* filehandle, uint64_t * RecNum,
                               uint64_t StopRecord, record_buf_t *buffer,
                               uint64_t *oflcorrection, uint64_t *timetag, int *channel)
{
    /*
     next_record() reads the next records of a file until it finds a photon, and then returns.
     Inputs:
     filehandle         FILE pointer with an open record file to read the photons
     RecNum             pointer to the index of the record being read
     StopRecord         Last record number of interest
     buffer             pointer to a record_buf_t structure which will be used for chunk file reading
     oflcorrection      pointer to an unsigned integer 64 bits. Will record the
                        time correction in the timetags due to overflow.
     timetag            pointer to an unsigned integer 64 bits. Timetag of the
                        next photon (see outputs for details).
     channel            pointer to an integer. Channel of the next photon (see outputs for details).
     Outputs:
     filehandle         FILE pointer with reader at the position of last analysed record
     RecNum             index of last analysed record
     oflcorrection      offset time on the timetags read in the file, due to overflows. Should not be used.
     timetag            timetag of the last photon read. It already includes the
                        overflow correction so the value can  be used directly.
     channel            channel of the last photon read. 0 will usually be
                        sync and >= 1 other input channels.
     Returns:
     1 when found a photon,
     0 when reached end of file.
     */
    pop_record:
    if (buffer->head < RECORD_CHUNK && (*RecNum)<StopRecord) { // still have records on buffer
        // This .c file is preprocessed by _readTTTRecords_build.py by
        // replacing the ##parser## tag with different parsers. This
        // replacing makes the file into a valid C file. By doing this
        // we can easily generate one library per record type. The ultimate
        // reason is to avoid using either a switch statment or calling a
        // a function via a function pointer inside a hot loop.
        Parse##parser##(buffer->records[buffer->head], channel, timetag, oflcorrection);
        buffer->head++;
        (*RecNum)++;

        return true;
    } else if (*RecNum >= StopRecord) { // run out of records
          return false;
    }
    // run out of buffer
    buffer->head = 0;
    load_buffer(buffer->records, filehandle);
    goto pop_record;
}



// = = = = = = = = = = =//
// TIME TRACE ALGORITHM //
// = = = = = = = = = = =//

typedef struct _timetrace_args {
        int end_of_header;
        int n_bins;
        int *ptr_trace;
        uint32_t *buffer;
        uint64_t *ptr_recnum;
        uint64_t RecNum_start;
        uint64_t RecNum_stop;
        int select_channel;
        uint64_t time_bin_length;
        char *filepath;
    } timetrace_args;

static inline THREAD_FUNC_DEF(timetrace_section) {
    /*Return an intensity time trace from a thread.
    Arguments:
    end_of_header:    Position of the end the file header in bytes
    n_bins:           Number of bins in a time trace
    ptr_trace:        Pointer to the memory where the timetrace is being stored
    ptr_recnum:       Pointer to the memory where the recnum trace is being stored
    RecNum_start:     First record number of interest
    RecNum_stop:      Last record number of interest
    time_bin_length:  Length of a time bin
    */
    // Get a filehandle local to the thread
    const timetrace_args *args = (timetrace_args*)arguments;

    // Open file and jump to first record
    FILE *filehandle = fopen(args->filepath, "rb");
    c_fseek(filehandle,
            (long int)(args->end_of_header + (4 * args->RecNum_start) ));

    // prepare record buffer
    record_buf_t TTTRRecord;
    TTTRRecord.records = args->buffer;
    record_buf_reset(&TTTRRecord);

    // return values for next photon
    uint64_t oflcorrection = 0;
    uint64_t timetag = 0;
    int channel = -1;

    bool record_arrived = true;
    uint64_t RecNum = args->RecNum_start;
    uint64_t end_of_bin;

    load_buffer(TTTRRecord.records, filehandle);
    int photon_counter=0;
    for (int i = 0; i < args->n_bins; i++)
    {
        end_of_bin = (i+1) * args->time_bin_length;
        photon_counter = 0;
        if (args->select_channel < 0) {
            while (timetag < end_of_bin && record_arrived) {
                record_arrived = next_record(filehandle, &RecNum, args->RecNum_stop,
                                             &TTTRRecord, &oflcorrection, &timetag, &channel);
                photon_counter += (channel >= 0);
            }
        } else {
            while (timetag < end_of_bin && record_arrived) {
                record_arrived = next_record(filehandle, &RecNum, args->RecNum_stop,
                                             &TTTRRecord, &oflcorrection, &timetag, &channel);
                photon_counter += (channel == args->select_channel);
            }
        }

        if (record_arrived) { // the last incomplete bin is discarded
            args->ptr_recnum[i] = RecNum;
            args->ptr_trace[i] = photon_counter;
        } else break; // no photons left
    }

    fclose(filehandle);
    return RETURN_SECTION;
}

void timetrace(char filepath[], int end_of_header, uint64_t RecNum_start,
               uint64_t NumRecords, uint64_t time_bin_length, int time_trace[],
               uint64_t RecNum_trace[], int select_channel, int nb_of_bins, int n_threads)
{
    int i, j, k; // looping indices

    // Make number or records a mutiple of the number of threads
    NumRecords = NumRecords-(NumRecords%n_threads);
    uint64_t records_per_thread = (uint64_t)(NumRecords/n_threads);

    // Prepare the threads
#if defined(__linux__ ) || defined(__APPLE__)
    pthread_t *tid_array;
    tid_array = (pthread_t*) malloc(n_threads * sizeof(pthread_t));
#elif defined(_WIN32)
    HANDLE *hThreadArray;
    hThreadArray = (HANDLE*) malloc(n_threads * sizeof(HANDLE));
    DWORD   *dwThreadIdArray;
    dwThreadIdArray = (DWORD*) malloc(n_threads * sizeof(DWORD));
#endif

    // Fill in the arguments for the different threads
    timetrace_args *thread_args;
    thread_args = (timetrace_args*) malloc(n_threads * sizeof(timetrace_args));

    for (i = 0; i < n_threads; ++i) {
        thread_args[i].ptr_trace = (int*) malloc(nb_of_bins * sizeof(int));
        thread_args[i].ptr_recnum = (uint64_t*) malloc(nb_of_bins * sizeof(uint64_t));
        for (j = 0; j < nb_of_bins; ++j) {
            // we are going to use the -1 as a flag
            // to find where we stopped adding bins
            thread_args[i].ptr_trace[j] = -1;
            thread_args[i].ptr_recnum[j] = -1;
        }

        thread_args[i].buffer = (uint32_t*) malloc(RECORD_CHUNK * sizeof(uint32_t));
        thread_args[i].end_of_header = end_of_header;
        thread_args[i].RecNum_start = (uint64_t)i * records_per_thread + RecNum_start;
        thread_args[i].RecNum_stop = ((uint64_t)i+1) * records_per_thread + RecNum_start;
        thread_args[i].n_bins = nb_of_bins;
        thread_args[i].select_channel = select_channel;
        thread_args[i].time_bin_length = time_bin_length;
        thread_args[i].filepath = filepath;
    }

    for (i = 0; i < n_threads; ++i) {
        CREATE_THREAD(timetrace_section);
    }

#if defined(__linux__ ) || defined(__APPLE__)
    for (int i = 0; i < n_threads; ++i) {
        pthread_join(tid_array[i], NULL);
    }
#elif defined(_WIN32)
    WaitForMultipleObjects(n_threads, hThreadArray, TRUE, INFINITE);
#endif


    // * = * = * = * = * = * = * = * = * = * = * = * = * = * =
    // Splice the timetraces of each thread into a single one.
    // Do also the recnums.
    // * = * = * = * = * = * = * = * = * = * = * = * = * = * =
    k = 0; // index for the time traces within a thread
    int new_val_tt;
    uint64_t new_val_rec; // stores the possible value for the trace
    j = 0;  // use j to go over the threads
    for (i = 0; i < nb_of_bins; ++i)
    {
        new_val_tt = thread_args[j].ptr_trace[k];
        new_val_rec = thread_args[j].ptr_recnum[k];

        if (new_val_tt != -1) {
            time_trace[i]=(int)new_val_tt;
            RecNum_trace[i]=new_val_rec;
            k++;
        } else { // found -1 time to go to next thread
            j++;
            if (j>=n_threads){ // we run out of threads
                break;
            }
            time_trace[i]=(int)thread_args[j].ptr_trace[0];
            RecNum_trace[i]=thread_args[j].ptr_recnum[0];
            k = 1;
        }
    }
    // Ended splicing partial timetraces

    // Free the memory and return
    for (i = 0; i < n_threads; ++i)
    {
        free(&(thread_args[i].ptr_trace[0]));
        free(&(thread_args[i].ptr_recnum[0]));
        free(&(thread_args[i].buffer[0]));
    }

#if defined(__linux__ ) || defined(__APPLE__)
    free(tid_array);
#elif defined(_WIN32)
    for (int i = 0; i < n_threads; ++i)
    {
        CloseHandle(hThreadArray[i]);
    }
    free(hThreadArray);
    free(dwThreadIdArray);
#endif
    free(thread_args);
    return;
}

// = = = = = = = //
// G2 ALGORITHMS //
// = = = = = = = //

enum mode {FAST, RING, CLASSIC, SYMMETRIC};

typedef struct _g2_args {
        int end_of_header;
        int n_bins;                // number of bins in histogram
        int first_range;
        int n_ranges;              // number of ranges in thread
        int channel_start;
        int channel_stop;
        size_t buffer_size;        // ring buffer size (only used by g2_ring)
        uint64_t correlation_window;  // duration of histogram
        int *ptr_hist;             // pointer to the thread histogram
        uint32_t *buffer;          // record buffer
        uint64_t *RecNum_start;    // array with the thread recnum starts
        uint64_t *RecNum_stop;     // array with the thread recnum stops
        char *filepath;
    } g2_args;

static inline THREAD_FUNC_DEF(g2_fast_section) {
    // Get a thread filehandle
    const g2_args *args = (g2_args*)arguments;

    // Prepare the file
    FILE *filehandle = fopen(args->filepath, "rb");

    record_buf_t TTTRRecord;
    TTTRRecord.records = args->buffer;
    record_buf_reset(&TTTRRecord);

    // return values next photon
    uint64_t oflcorrection = 0;
    uint64_t start_time;
    uint64_t stop_time;
    int channel = -1;

    // variables for g2 algo
    const int channel_start = args->channel_start;
    const int channel_stop = args->channel_stop;

    uint64_t i = 0;
    uint64_t delta=0;
    uint64_t correlation_window = args->correlation_window;
    const int nb_of_bins = args->n_bins;
    const uint64_t resolution = correlation_window / nb_of_bins;

    uint64_t RecNum;
    uint64_t RecNum_STOP;
    bool record_arrived = true;

    // loop over postselection ranges assigned to thread
    for (int range_idx = 0; range_idx < args->n_ranges; range_idx++) {
        record_arrived = true;
        RecNum = args->RecNum_start[args->first_range + range_idx];
        RecNum_STOP = args->RecNum_stop[args->first_range + range_idx];
        c_fseek(filehandle,
               (long int)(args->end_of_header + (4 * RecNum)));

        // start g2 algo
        // prefill record buffer
        load_buffer(TTTRRecord.records, filehandle);
        TTTRRecord.head = 0;
        channel = -1;

        while(record_arrived){
            // FIND NEXT START PHOTON
            while(record_arrived && channel != channel_start){
                record_arrived = next_record(filehandle, &RecNum, RecNum_STOP,
                                             &TTTRRecord, &oflcorrection, &start_time, &channel);
            }
            // found start photon

            // FIND NEXT STOP PHOTON
            while (record_arrived && channel != channel_stop) {
                record_arrived = next_record(filehandle, &RecNum, RecNum_STOP,
                                             &TTTRRecord, &oflcorrection, &stop_time, &channel);
            }
            // found stop photon

            // ADD DELAY TO HISTOGRAM
            delta = stop_time - start_time;
            if (delta < correlation_window && record_arrived) {
                i = delta / resolution;
                args->ptr_hist[i]++;
            }
        } // end g2 algo
    }
    fclose(filehandle);
    return RETURN_SECTION;
}

static inline THREAD_FUNC_DEF(g2_symmetric_section) {
    // Get a thread filehandle
    const g2_args *args = (g2_args*)arguments;

    // Prepare the file
    FILE *filehandle = fopen(args->filepath, "rb");

    // prepare record buffer
    record_buf_t TTTRRecord;
    TTTRRecord.records = args->buffer;
    record_buf_reset(&TTTRRecord);

    // return values next photon
    uint64_t oflcorrection = 0;
    int channel = -1;
    uint64_t timetag = 0;

    // variables for g2 algo
    uint64_t delta, idx;
    int *ptr_hist = args->ptr_hist;
    const uint64_t nb_of_bins = args->n_bins;
    const uint64_t central_bin = nb_of_bins/2;
    const int channel_start = args->channel_start;
    const int channel_stop = args->channel_stop;
    uint64_t correlation_window = args->correlation_window;
    const uint64_t resolution = correlation_window / nb_of_bins;

    int i;  // index for the loop over ring buffer
    uint64_t RecNum;
    uint64_t RecNum_STOP;

    // Prepare the ring buffer for the start photons
    ring_buf_t cbuf_2 = ring_buf_allocate((int)args->buffer_size);
    ring_buf_t cbuf_1 = ring_buf_allocate((int)args->buffer_size);

    // loop over the postselection ranges assigned to thread
    for (int range_idx = 0; range_idx < args->n_ranges; range_idx++) {
        RecNum = args->RecNum_start[args->first_range + range_idx];
        RecNum_STOP = args->RecNum_stop[args->first_range + range_idx];
        c_fseek(filehandle,
            (long int)(args->end_of_header + (4 * RecNum) ));

        // start g2 algo
        // prefill ring buffer
        load_buffer(TTTRRecord.records, filehandle);
        TTTRRecord.head = 0;
        while(next_record(filehandle, &RecNum, RecNum_STOP, &TTTRRecord,
                          &oflcorrection, &timetag, &channel)) {

            if (channel == channel_start) {
                ring_buf_put(&cbuf_1, timetag);
                check_and_grow_buf(&cbuf_1, timetag, correlation_window);
                for(i = cbuf_2.head-1; i > (cbuf_2.head-1-cbuf_2.count); i--) {
                    delta = timetag - cbuf_2.buffer[(i+2*cbuf_2.count)%cbuf_2.count];
                    idx = central_bin - delta / resolution - 1;
                    if (delta < correlation_window && idx < nb_of_bins) {
                        ptr_hist[idx]++;
                    } else break;
                }
                continue;
            }

            if (channel == channel_stop) {
                ring_buf_put(&cbuf_2, timetag);
                check_and_grow_buf(&cbuf_2, timetag, correlation_window);
                for(i = cbuf_1.head-1; i > (cbuf_1.head-1-cbuf_1.count); i--) {
                    delta = timetag - cbuf_1.buffer[(i+2*cbuf_1.count)%cbuf_1.count];
                    idx = central_bin + delta / resolution;
                    if (delta < correlation_window && idx < nb_of_bins) {
                        ptr_hist[idx]++;
                    } else break;
                }
            }
        } // end g2 algo
    }
    free(cbuf_1.buffer);
    free(cbuf_2.buffer);
    fclose(filehandle);
    return RETURN_SECTION;
}

static inline THREAD_FUNC_DEF(g2_ring_section) {
    // Get a thread filehandle
    const g2_args *args = (g2_args*)arguments;

    // Prepare the file
    FILE *filehandle = fopen(args->filepath, "rb");

    // prepare record buffer
    record_buf_t TTTRRecord;
    TTTRRecord.records = args->buffer;
    record_buf_reset(&TTTRRecord);

    // return values next photon
    uint64_t oflcorrection = 0;
    int channel = -1;
    uint64_t timetag = 0;

    // variables for g2 algo
    uint64_t delta, idx;
    int *ptr_hist = args->ptr_hist;
    const int nb_of_bins = args->n_bins;
    const int channel_start = args->channel_start;
    const int channel_stop = args->channel_stop;
    uint64_t correlation_window = args->correlation_window;
    const uint64_t resolution = correlation_window / nb_of_bins;

    int i;  // index for the loop over ring buffer
    uint64_t RecNum;
    uint64_t RecNum_STOP;

    // Prepare the ring buffer for the start photons
    ring_buf_t cbuf = ring_buf_allocate((int)args->buffer_size);

    // loop over the postselection ranges assigned to thread
    for (int range_idx = 0; range_idx < args->n_ranges; range_idx++) {
        ring_buf_reset(&cbuf);
        RecNum = args->RecNum_start[args->first_range + range_idx];
        RecNum_STOP = args->RecNum_stop[args->first_range + range_idx];
        c_fseek(filehandle,
            (long int)(args->end_of_header + (4 * RecNum) ));

        // start g2 algo
        // prefill ring buffer
        load_buffer(TTTRRecord.records, filehandle);
        TTTRRecord.head = 0;
        while(next_record(filehandle, &RecNum, RecNum_STOP, &TTTRRecord,
                          &oflcorrection, &timetag, &channel)) {

            if (channel == channel_start) {
                ring_buf_put(&cbuf, timetag);
                check_and_grow_buf(&cbuf, timetag, correlation_window);
                continue;
            }

            if (channel == channel_stop && cbuf.count > 0) {
                for(i = cbuf.head-1; i > (cbuf.head-1-cbuf.count); i--) {
                    delta = timetag - cbuf.buffer[(i+cbuf.count)%cbuf.count];
                    if (delta < correlation_window) {
                        idx = delta / resolution;
                        ptr_hist[idx]++;
                    } else break;
                }
            }
        } // end g2 algo
    }
    free(cbuf.buffer);
    fclose(filehandle);
    return RETURN_SECTION;
}

static inline THREAD_FUNC_DEF(g2_classic_section) {
    // Get a thread filehandle
    const g2_args *args = (g2_args*)arguments;

    // Prepare the file
    FILE *filehandle = fopen(args->filepath, "rb");

    record_buf_t TTTRRecord;
    TTTRRecord.records = args->buffer;
    record_buf_reset(&TTTRRecord);

    node_t* start_buff_head = NULL;
    node_t* stop_buff_head = NULL;
    int start_buff_length = 0;
    int stop_buff_length = 0;
    node_t* stop_corr_buff_head = NULL;
    int stop_corr_buff_length = 0;
    node_t* current = NULL;
    uint64_t correlation_window_end = 0;
    uint64_t start_time = 0;
    uint64_t oflcorrection = 0;
    uint64_t timetag = 0;
    int channel = -1;
    uint64_t i = 0;
    const uint64_t correlation_window = args->correlation_window;

    const int channel_start = args->channel_start;
    const int channel_stop = args->channel_stop;
    uint64_t RecNum, RecNum_STOP;
    const int nb_of_bins = args->n_bins;
    bool record_arrived;

    for (int range_idx = 0; range_idx < args->n_ranges; range_idx++) {
        // reset file reader and go to the start position RecNum_start
        record_arrived=1;
        RecNum = args->RecNum_start[args->first_range + range_idx];
        RecNum_STOP = args->RecNum_stop[args->first_range + range_idx];
        c_fseek(filehandle,
            (long int)(args->end_of_header + (4 * RecNum) ));
        // First item in the chained lists will be kept as anchor and only the 'next' items will contain timetags.
        // This avoids emptying totally the list and having to recreate it when starting to fill it again.
        head_init(&start_buff_head, &start_buff_length);
        head_init(&stop_buff_head, &stop_buff_length);
        head_init(&stop_corr_buff_head, &stop_corr_buff_length);

        /*
         This algorithm implies using 3 buffers:
         start_buff_head      : the start photons buffer, where all unused start photons go (to be used later)
         stop_buff_head       : the stop photons buffer, where all unused stop photons go (to be used later)
         stop_corr_buff_head  : the correlation stop photons buffer. This buffer contains all the stop photons which fit
         in a correlation window around the selected start photon. For each new start photon,
         it needs to be modified removing the old photons which do not fit anymore in the correlation
         window and adding the new ones which now fit in the correlation window.

         Note that this algorithm supposes the list of photons to be ordered chronologically.
         */

        // while there are still unread photons in the file or unused start photons in the buffer
        // prefill the TTTRRecord struct
        load_buffer(TTTRRecord.records, filehandle);
        while(record_arrived || start_buff_length > 0){
            // FIND NEXT START PHOTON
            // first, take first start photon in buffer
            if(start_buff_length > 0){
                start_time = pop(start_buff_head, &start_buff_length);
            }
            // if start buffer is empty, read photons until a start photon is found, and feed stop buffer in the process
            else {
                channel = -1;
                while(channel != channel_start && record_arrived){
                    record_arrived = next_record(filehandle, &RecNum, RecNum_STOP,
                              &TTTRRecord, &oflcorrection, &timetag, &channel);
                    if (channel == channel_stop){ // store in stop photons buffer
                        push(stop_buff_head, timetag, &stop_buff_length);
                    }
                    else if (channel == channel_start) { // channel 0
                        start_time = timetag;
                    }
                }
                if (channel != channel_start && record_arrived) {
                    break;
                }
            }
            correlation_window_end = start_time + correlation_window;

            // FIND ALL STOP PHOTONS IN CORRELATION WINDOW
            // complete stop photons array with new stop photons from buffer fitting in correlation window
            while(stop_buff_length > 0 && stop_buff_head->next->val < correlation_window_end) {
                push(stop_corr_buff_head, pop(stop_buff_head, &stop_buff_length), &stop_corr_buff_length);
            }

            // if stop buffer is empty, read photons until the time gets out of the
            // correlation window, and feed start buffer and the stop photons array in the process
            if (stop_buff_length == 0) {
                while (timetag < correlation_window_end && record_arrived) {
                    record_arrived = next_record(filehandle, &RecNum, RecNum_STOP,
                              &TTTRRecord, &oflcorrection, &timetag, &channel);
                    // start photon -> store in start photon buffer (to be used later)
                    if (channel == channel_start) {
                        push(start_buff_head, timetag, &start_buff_length);
                    }
                    // stop photon
                    else if (channel == channel_stop) {
                        // qualifies in correlation window -> store in correlation window stop buffer
                        if (timetag < correlation_window_end) {
                            push(stop_corr_buff_head, timetag, &stop_corr_buff_length);
                        }
                        // doesn't qualify -> store in stop photon buffer (to be used later)
                        else {
                            push(stop_buff_head, timetag, &stop_buff_length);
                        }
                    }
                }
            }
            // remove stop photons which are out of the correlation window (pop)
            for(i = 0; i < (uint64_t) stop_corr_buff_length; i++) {
                if(stop_corr_buff_head->next->val < start_time) {
                    pop(stop_corr_buff_head, &stop_corr_buff_length);
                }
                else {
                    break;
                }
            }
            // perform a histogram of the stop times - start time and add it to the main histogram result
            current = stop_corr_buff_head->next;
            while(current != NULL && record_arrived) {
                if (current->val - start_time < correlation_window) {
                    i = (uint64_t) (current->val - start_time) * nb_of_bins / correlation_window;
                    args->ptr_hist[i]++;
                }
                current = current->next;
            }
        }
    // When we are done we have to clear the memory for the linked list
    while(pop(start_buff_head, &start_buff_length)){}
    free(start_buff_head);
    while(pop(stop_buff_head, &stop_buff_length)){}
    free(stop_buff_head);
    while(pop(stop_corr_buff_head, &stop_corr_buff_length)){}
    free(stop_corr_buff_head);
    }

    fclose(filehandle);
    return RETURN_SECTION;
}

void calculate_g2(char filepath[], int end_of_header,
                  uint64_t *RecNum_start, uint64_t *RecNum_stop,
                  int nb_of_ranges, uint64_t max_time, int histogram[],
                  int nb_of_bins, int channel_start, int channel_stop,
                  int buffer_size, int n_threads, int mode)
{
    // Variables used in algo to distribute ranges evenly among threads
    int first_range = 0;
    int n_ranges;
    int leftover_ranges = nb_of_ranges;

    if (nb_of_ranges < n_threads){
        n_threads = nb_of_ranges;

        printf("%s\n", "There are less post selection ranges than threads.");
        printf("%s\n", "Computation will continue running with one thread per range.");
    }

    // Prepare the threads
#if defined(__linux__ ) || defined(__APPLE__)
    pthread_t *tid_array;
    tid_array = (pthread_t*) malloc(n_threads * sizeof(pthread_t));
#elif defined(_WIN32)
    HANDLE *hThreadArray;
    hThreadArray = (HANDLE*) malloc(n_threads * sizeof(HANDLE));
    DWORD   *dwThreadIdArray;
    dwThreadIdArray = (DWORD*) malloc(n_threads * sizeof(DWORD));
#endif

    // Fill in the arguments for the different threads
    g2_args *thread_args;
    thread_args = (g2_args*) malloc(n_threads * sizeof(g2_args));

    for (int i = 0; i < n_threads; ++i) {
        thread_args[i].ptr_hist= (int*) malloc(nb_of_bins * sizeof(int));
        for (int j = 0; j < nb_of_bins; ++j) {
            thread_args[i].ptr_hist[j] = 0;
        }

        thread_args[i].end_of_header = end_of_header;
        thread_args[i].n_bins = nb_of_bins;
        thread_args[i].correlation_window = max_time;

        thread_args[i].buffer = (uint32_t*) malloc(RECORD_CHUNK * sizeof(uint32_t));
        thread_args[i].buffer_size = buffer_size;

        thread_args[i].RecNum_start = RecNum_start;
        thread_args[i].RecNum_stop = RecNum_stop;

        thread_args[i].channel_start = channel_start;
        thread_args[i].channel_stop = channel_stop;

        // Simple algorithm to distribute ranges among threads
        thread_args[i].first_range = first_range;
        n_ranges = (leftover_ranges / (n_threads-i)) +
                   (leftover_ranges % (n_threads-i) ? 1 : 0);
        thread_args[i].n_ranges = n_ranges;
        leftover_ranges -= n_ranges;
        first_range += n_ranges;


        thread_args[i].filepath = filepath;
    }

    // Print what record number ranges are associated with each thread
    for (int i = 0; i < n_threads; ++i)
    {
        printf("Thread:%d\n", i );
        for(int j = 0; j < thread_args[i].n_ranges; j++) {
            printf("[%llu, %llu]\n", thread_args[i].RecNum_start[thread_args[i].first_range + j],
                                      thread_args[i].RecNum_stop[thread_args[i].first_range + j]);
        }
    }

    printf("%s%d\n", "G2 mode: ", mode);
    for (int i = 0; i < n_threads; ++i) {
        switch (mode) {
            case FAST:
                CREATE_THREAD(g2_fast_section);
                break;
            case RING:
                CREATE_THREAD(g2_ring_section);
                break;
            case CLASSIC:
                CREATE_THREAD(g2_classic_section);
                break;
            case SYMMETRIC:
                CREATE_THREAD(g2_symmetric_section);
                break;
            default:
                printf("%s\n", "NON-EXISTENT G2 MODE");
                goto free_memory;
        }
    }

#if defined(__linux__ ) || defined(__APPLE__)
    for (int i = 0; i < n_threads; ++i) {
        pthread_join(tid_array[i], NULL);
    }
#elif defined(_WIN32)
    WaitForMultipleObjects(n_threads, hThreadArray, TRUE, INFINITE);
#endif

    // Combine the histograms
    for (int i = 0; i < n_threads; ++i)
    {
        for (int j = 0; j < nb_of_bins; ++j)
        {
            histogram[j] += thread_args[i].ptr_hist[j];
        }
    }

    free_memory:
    for (int i = 0; i < n_threads; ++i)
    {
        free(&(thread_args[i].ptr_hist[0]));
        free(&(thread_args[i].buffer[0]));
    }

#if defined(__linux__ ) || defined(__APPLE__)
    free(tid_array);
#elif defined(_WIN32)
    for (int i = 0; i < n_threads; ++i)
    {
        CloseHandle(hThreadArray[i]);
    }
    free(hThreadArray);
    free(dwThreadIdArray);
#endif
    free(thread_args);
    return;
}