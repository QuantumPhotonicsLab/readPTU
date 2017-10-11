//
//  main.c
//  ptu-read-tests
//
//  Created by Raphaël Proux on 03/10/2017.
//  Copyright © 2017 Raphaël Proux. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>

// RecordTypes
#define rtPicoHarpT3     0x00010303    // (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $03 (PicoHarp)
#define rtPicoHarpT2     0x00010203    // (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $03 (PicoHarp)
#define rtHydraHarpT3    0x00010304    // (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $03 (T3), HW: $04 (HydraHarp)
#define rtHydraHarpT2    0x00010204    // (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $04 (HydraHarp)
#define rtHydraHarp2T3   0x01010304    // (SubID = $01 ,RecFmt: $01) (V2), T-Mode: $03 (T3), HW: $04 (HydraHarp)
#define rtHydraHarp2T2   0x01010204    // (SubID = $01 ,RecFmt: $01) (V2), T-Mode: $02 (T2), HW: $04 (HydraHarp)
#define rtTimeHarp260NT3 0x00010305    // (SubID = $00 ,RecFmt: $01) (V2), T-Mode: $03 (T3), HW: $05 (TimeHarp260N)
#define rtTimeHarp260NT2 0x00010205    // (SubID = $00 ,RecFmt: $01) (V2), T-Mode: $02 (T2), HW: $05 (TimeHarp260N)
#define rtTimeHarp260PT3 0x00010306    // (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T3), HW: $06 (TimeHarp260P)
#define rtTimeHarp260PT2 0x00010206    // (SubID = $00 ,RecFmt: $01) (V1), T-Mode: $02 (T2), HW: $06 (TimeHarp260P)

// this 'node' structure is a chained list to be used with the buffers. It is easy to add an item at the end and read+remove an item at the beginning (see the three functions below).
typedef struct node {
    uint64_t val;
    struct node * next;
} node_t;

void head_init(node_t** head, int* length) {
    // initialize a chained list with value 0 and length 0.
    *head = malloc(sizeof(node_t));
    (*head)->val = 0;
    (*head)->next = NULL;
    *length = 0;
}

void push(node_t * head, uint64_t val, int* length) {
    // add an item at the end of the list
    node_t * current = head;
    while (current->next != NULL) {
        current = current->next;
    }
    
    // now we can add a new variable
    current->next = malloc(sizeof(node_t));
    current->next->val = val;
    current->next->next = NULL;
    *length = *length + 1;
}

uint64_t pop(node_t * head, int* length) {
    // remove the first info item (second in the list) from the list, returning its value
    // remember that, for simplicity, we keep the very first item of the list alive so we don't need to recreate one each time the list becomes empty.
    uint64_t retval = 0;
    node_t * next_node = NULL;
    
    if (head->next == NULL) {
        return 0;
    }
    
    next_node = head->next->next;
    retval = head->next->val;
    free(head->next);
    head->next = next_node;
    *length = *length - 1;
    
    return retval;
}


int c_fseek(FILE *filehandle, long int offset)
{
    return fseek(filehandle, offset, SEEK_SET);
}


void ProcessPHT2(FILE* filehandle, uint64_t *oflcorrection, uint64_t *timetag, int *channel)
{
    /*
     ProcessPHT2() reads the next records of a file until it finds a photon, and then returns.
     Inputs:
     filehandle         FILE pointer with an open record file to read the photons
     oflcorrection      pointer to an unsigned integer 64 bits. Will record the time correction in the timetags due to overflow (see output for details).
     timetag            pointer to an unsigned integer 64 bits (see outputs for details).
     channel            pointer to an integer (see outputs for details).
     Outputs:
     filehandle         FILE pointer with reader at the position of last analysed record
     oflcorrection      offset time on the timetags read in the file, due to overflows.
     timetag            if a photon is read, timetag of this photon. Otherwise, timetag == 0. It already includes the overflow correction so the value can be used directly.
     channel            if a photon is read, channel of this photon. 0 will usually be sync and >= 1 other input channels. If the record is not a photon, channel == -1 for an overflow record, -2 for a marker record.
     */
    /* SHOULD WORK BUT FUNCTION NOT TESTED */
    
    const int T2WRAPAROUND = 210698240;
    union
    {
        unsigned int allbits;
        struct
        {
            unsigned time   :28;
            unsigned channel  :4;
        } bits;
        
    } Record;
    unsigned int markers;
    uint32_t TTTRRecord = 0;
    
    fread(&TTTRRecord, 1, sizeof(TTTRRecord) ,filehandle);
    Record.allbits = TTTRRecord;
    
    if(Record.bits.channel == 0xF) //this means we have a special record
    {
        //in a special record the lower 4 bits of time are marker bits
        markers = Record.bits.time & 0xF;
        if(markers == 0) //this means we have an overflow record
        {
            *oflcorrection += T2WRAPAROUND; // unwrap the time tag overflow
            *channel = -1;
        }
        else //a marker
        {
            //Strictly, in case of a marker, the lower 4 bits of time are invalid
            //because they carry the marker bits. So one could zero them out.
            //However, the marker resolution is only a few tens of nanoseconds anyway,
            //so we can just ignore the few picoseconds of error.
            *timetag = *oflcorrection + Record.bits.time;
            *channel = -2;
        }
    }
    else
    {
        if((int)Record.bits.channel > 4) //Should not occur
        {
            *timetag = 0;
            *channel = -3;
        }
        else
        {
            *timetag = *oflcorrection + Record.bits.time;
            *channel = Record.bits.channel;
        }
    }
}

void ProcessHHT2(FILE* filehandle, int HHVersion, uint64_t *oflcorrection, uint64_t *timetag, int *channel)
{
    /*
     ProcessHHT2() reads the next records of a file until it finds a photon, and then returns.
     Inputs:
     filehandle         FILE pointer with an open record file to read the photons
     HHVersion          Hydrahard version 1 or 2. Depends on record type specification in the header.
     oflcorrection      pointer to an unsigned integer 64 bits. Will record the time correction in the timetags due to overflow (see output for details).
     timetag            pointer to an unsigned integer 64 bits (see outputs for details).
     channel            pointer to an integer (see outputs for details).
     Outputs:
     filehandle         FILE pointer with reader at the position of last analysed record
     oflcorrection      offset time on the timetags read in the file, due to overflows.
     timetag            if a photon is read, timetag of this photon. Otherwise, timetag == 0. It already includes the overflow correction so the value can be used directly.
     channel            if a photon is read, channel of this photon. 0 will usually be sync and >= 1 other input channels. If the record is not a photon, channel == -1 for an overflow record, -2 for a marker record.
     */
    /* FUNCTION TESTED */
    const uint64_t T2WRAPAROUND_V1 = 33552000;
    const uint64_t T2WRAPAROUND_V2 = 33554432;
    union{
        uint32_t   allbits;
        struct{ unsigned timetag  :25;
            unsigned channel  :6;
            unsigned special  :1; // or sync, if channel==0
        } bits;
    } T2Rec;
    uint32_t TTTRRecord = 0;
    
    
    fread(&TTTRRecord, 1, sizeof(TTTRRecord) ,filehandle);
    
    T2Rec.allbits = TTTRRecord;
    
    if(T2Rec.bits.special==1)
    {
        if(T2Rec.bits.channel==0x3F) //an overflow record
        {
            if(HHVersion == 1)
            {
                *oflcorrection += T2WRAPAROUND_V1;
                *timetag = 0;
                *channel = -1;
            }
            else
            {
                //number of overflows is stored in timetag
                if(T2Rec.bits.timetag==0) //if it is zero it is an old style single overflow
                {
                    *oflcorrection += T2WRAPAROUND_V2;  //should never happen with new Firmware!
                }
                else
                {
                    *oflcorrection += T2WRAPAROUND_V2 * T2Rec.bits.timetag;
                }
            }
            
            *timetag = 0;
            *channel = -1;
        }
        
        if((T2Rec.bits.channel>=1)&&(T2Rec.bits.channel<=15)) //markers
        {
            *timetag = *oflcorrection + T2Rec.bits.timetag;
            //Note that actual marker tagging accuracy is only some ns.
            *channel = -2;
            // *channel = T2Rec.bits.channel;
        }
        
        else if(T2Rec.bits.channel==0) //sync
        {
            *timetag = *oflcorrection + T2Rec.bits.timetag;
            *channel = T2Rec.bits.channel;
        }
    }
    else //regular input channel
    {
        *timetag = *oflcorrection + T2Rec.bits.timetag;
        *channel = T2Rec.bits.channel + 1;
    }
    //    printf("%d %llu\n", *channel, *timetag);
}


void RecordHHT2(FILE* filehandle)
{
    /*
     RecordHHT2() is a function written for test purposes only. It writes a dummy file containing prearranged photons to test the timetrace and g2 algorithms.
     */
    union {
        uint32_t   allbits;
        struct{ unsigned timetag  :25;
            unsigned channel  :6;
            unsigned special  :1; // or sync, if channel==0
        } bits;
    } T2Rec;
    uint32_t TTTRRecord = 0;
    int i = 0;
    for(i = 0; i < 100; i++){
        if (i % 2){
            T2Rec.bits.special = 0;  // regular input channel
            T2Rec.bits.channel = 0;  // channel number
            T2Rec.bits.timetag = i * 500 + 10;  // timetag in ps
        }
        else {
            T2Rec.bits.special = 1;  // sync channel
            T2Rec.bits.channel = 0;  // sync channel
            T2Rec.bits.timetag = i * 500 + 10;  // timetag in ps
        }
        TTTRRecord = T2Rec.allbits;
        //        printf("%u\n",TTTRRecord);
        fwrite(&TTTRRecord, 4, 1, filehandle);
    }
}


int next_photon(FILE* filehandle, long long record_type, uint64_t *RecNum, uint64_t NumRecords, uint64_t *oflcorrection, uint64_t *timetag, int *channel)
{
    /*
     next_photon() reads the next records of a file until it finds a photon, and then returns.
     Inputs:
     filehandle         FILE pointer with an open record file to read the photons
     record_type        record type which depends on the device which recorded the file (see constants at the beginning of file)
     RecNum             pointer to the index of the record being read
     NumRecords         total number of records
     oflcorrection      pointer to an unsigned integer 64 bits. Will record the time correction in the timetags due to overflow. At the start of the file should be provided as 0 for initial value.
     timetag            pointer to an unsigned integer 64 bits. Timetag of the next photon (see outputs for details).
     channel            pointer to an integer. Channel of the next photon (see outputs for details).
     Outputs:
     filehandle         FILE pointer with reader at the position of last analysed record
     RecNum             index of last analysed record
     oflcorrection      offset time on the timetags read in the file, due to overflows. Should not be used.
     timetag            timetag of the last photon read. It already includes the overflow correction so the value can  be used directly.
     channel            channel of the last photon read. 0 will usually be sync and >= 1 other input channels.
     Returns:
     1 when found a photon,
     0 when reached end of file.
     */
    do {
        if(*RecNum < NumRecords) {
            switch (record_type) {
                case rtPicoHarpT2:
                    ProcessPHT2(filehandle, oflcorrection, timetag, channel);
                    *RecNum = *RecNum + 1;
                    break;
                case rtPicoHarpT3:
                    //ProcessPHT3(TTTRRecord);
                    //*RecNum = *RecNum + 1;
                    break;
                case rtHydraHarpT2:
                    ProcessHHT2(filehandle, 1, oflcorrection, timetag, channel);
                    *RecNum = *RecNum + 1;
                    break;
                case rtHydraHarpT3:
                    //ProcessHHT3(TTTRRecord, 1);
                    //*RecNum = *RecNum + 1;
                    break;
                case rtHydraHarp2T2:
                case rtTimeHarp260NT2:
                case rtTimeHarp260PT2:
                    ProcessHHT2(filehandle, 2, oflcorrection, timetag, channel);
                    *RecNum = *RecNum + 1;
                    break;
                case rtHydraHarp2T3:
                case rtTimeHarp260NT3:
                case rtTimeHarp260PT3:
                    //ProcessHHT3(TTTRRecord, 2);
                    //*RecNum = *RecNum + 1;
                    break;
                default:
                    break;
            }
        }
        else {
            break;
        }
    } while (*channel < 0); // while not a photon
    if (*RecNum >= NumRecords) {
        return 0;  // reached end of records (end of file)
    }
    else {
        return 1;  // found a photon
    }
}

void timetrace(FILE* filehandle, long long record_type, int end_of_header, uint64_t *RecNum, uint64_t NumRecords, uint64_t time_bin_length, uint64_t *time_vector, int *time_trace, uint64_t *RecNum_trace, int nb_of_bins)
{
    /*
     timetrace() computes the timetrace of a given measurement file. It does not differentiate the channel detecting the photons, which are all added together.
     Inputs:
     filehandle         FILE pointer with an open record file to read the photons
     record_type        record type which depends on the device which recorded the file (see constants at the beginning of file)
     end_of_header      offset in bytes to the beginning of the record section in the file
     RecNum             pointer to the index of the record being read
     NumRecords         total number of records
     time_bin_length    length of a time bin of the timetrace in picoseconds
     time_vector        preallocated array used for the x-axis of the timetrace. Should have nb_of_bins elements.
     time_trace         preallocated array used for the timetrace value. Should have nb_of_bins elements.
     RecNum_trace       preallocated array used for the values of RecNum corresponding to the last photon of each time bin. Should have nb_of_bins elements.
     nb_of_bins         number of bins for the timetrace (should correspond to the length of the time_trace array)
     Outputs:
     filehandle         FILE pointer with reader at the position of last analysed record
     RecNum             index of last analysed record
     time_trace         calculated timetrace
     */
    // IMPORTANT NOTE: every time in picoseconds
    uint64_t oflcorrection = 0;
    uint64_t timetag = 0;
    int channel = -1;
    uint64_t end_of_bin = 0;
    int add_photon_to_next_bin = 0;
    int i = 0;
    
    // reset file reader
    c_fseek(filehandle, end_of_header);
    
    for (i = 0; i < nb_of_bins; i++)
    {
        time_vector[i] = i * time_bin_length;
        end_of_bin = time_vector[i] + time_bin_length;
        // the last timetag read is in the bin
        if (timetag < end_of_bin) {
            // if we are starting (still didn't read a photon), add_photon_to_next_bin == 0, otherwise, == 1
            time_trace[i] = add_photon_to_next_bin;
            add_photon_to_next_bin = 0;
            while(*RecNum < NumRecords)
            {
                next_photon(filehandle, record_type, RecNum, NumRecords, &oflcorrection, &timetag, &channel);
                // if this photon is in the current bin
                if(timetag < end_of_bin)
                {
                    time_trace[i] = time_trace[i] + 1;
                }
                // else in belongs to some further bin (CAUTION: may not be the one immediately after)
                else
                {
                    RecNum_trace[i] = *RecNum;
                    add_photon_to_next_bin = 1;
                    break;
                }
            }
            if (*RecNum >= NumRecords) {  // for the last time bin
                RecNum_trace[i] = *RecNum;
            }
        }
        // not photon found in this time bin.
        else {
            time_trace[i] = 0;
        }
    }
}

void calculate_g2(FILE* filehandle, long long record_type, int end_of_header, uint64_t *RecNum, uint64_t NumRecords, uint64_t RecNum_start, uint64_t RecNum_stop, uint64_t *time_vector, int *histogram, int nb_of_bins, int channel_start, int channel_stop)
{
    /*
     calculate_g2() computes the g2 directly reading the measurement file. It uses a more complex algorithm than calculate_g2_fast(). This function will keep all photons in memory buffers, such that each start photon will be measured in regard of all the stop photons detected in a correlation window around it. This way, the measurement does not stop at the first stop photon but will take into account longer time scales. It is therefore safer to use with high photon count rates.
     Inputs:
     filehandle         FILE pointer with an open record file to read the photons
     record_type        record type which depends on the device which recorded the file (see constants at the beginning of file)
     end_of_header      offset in bytes to the beginning of the record section in the file
     RecNum             pointer to the index of the record being read
     NumRecords         total number of records
     RecNum_start       start of the section of records to analyse for the g2 (in terms of record index)
     RecNum_stop        stop of the section of records to analyse for the g2
     time_vector        precalculated array of times used for the x-axis of the g2 histogram. Should have nb_of_bins + 1 elements.
     histogram          preallocated array of zeros used for the g2 histogram. Should have nb_of_bins elements.
     nb_of_bins         number of bins for the histogram (should correspond to the length of the histogram array)
     channel_start      channel number used for start photons (sync will generally be 0)
     channel_stop       channel number used for stop photons (> 0, often 1)
     Outputs:
     filehandle         FILE pointer with reader at the position of last analysed record
     RecNum             index of last analysed record
     histogram          calculated g2 histogram
     */
    
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
    int i = 0;
    uint64_t correlation_window = 0;
    //    long next_print = 0;
    correlation_window = time_vector[nb_of_bins];
    
    // First item in the chained lists will be kept as anchor and only the 'next' items will contain timetags.
    // This avoids emptying totally the list and having to recreate it when starting to fill it again.
    head_init(&start_buff_head, &start_buff_length);
    head_init(&stop_buff_head, &stop_buff_length);
    head_init(&stop_corr_buff_head, &stop_corr_buff_length);
    
    /*
     This algorithm implies using 3 buffers:
     start_buff_head      : the start photons buffer, where all unused start photons go (to be used later)
     stop_buff_head       : the stop photons buffer, where all unused stop photons go (to be used later)
     stop_corr_buff_head  : the correlation stop photons buffer. This buffer contains all the stop photons which fit in a correlation window around the selected start photon. For each new start photon, it needs to be modified removing the old photons which do not fit anymore in the correlation window and adding the new ones which now fit in the correlation window.
     
     Note that this algorithm supposes the list of photons to be ordered chronologically.
     */
    
    // reset file reader and go to the start position RecNum_start
    c_fseek(filehandle, end_of_header + 4 * RecNum_start);
    *RecNum = RecNum_start;
    
    // while there are still unread photons in the file or unused start photons in the buffer
    while((*RecNum < RecNum_stop && *RecNum < NumRecords) || start_buff_length > 0){
        //        if (*RecNum > next_print){
        //            printf("%ld/%ld\n", *RecNum, RecNum_stop);
        //            next_print = next_print + 1000000;
        //        }
        
        // FIND NEXT START PHOTON
        // first, take first start photon in buffer
        if(start_buff_length > 0){
            start_time = pop(start_buff_head, &start_buff_length);
        }
        // if start buffer is empty, read photons until a start photon is found, and feed stop buffer in the process
        else {
            channel = -1;
            while(channel != channel_start && (*RecNum < RecNum_stop && *RecNum < NumRecords)){
                next_photon(filehandle, record_type, RecNum, NumRecords, &oflcorrection, &timetag, &channel);
                if (channel == channel_stop){ // store in stop photons buffer
                    push(stop_buff_head, timetag, &stop_buff_length);
                }
                else { // channel 0
                    start_time = timetag;
                }
            }
            if (channel != channel_start && (*RecNum >= RecNum_stop || *RecNum >= NumRecords)) {
                break;
            }
        }
        correlation_window_end = start_time + correlation_window;
        
        // FIND ALL STOP PHOTONS IN CORRELATION WINDOW
        // complete stop photons array with new stop photons from buffer fitting in correlation window
        while(stop_buff_length > 0 && stop_buff_head->next->val < correlation_window_end) {
            push(stop_corr_buff_head, pop(stop_buff_head, &stop_buff_length), &stop_corr_buff_length);
        }
        
        // if stop buffer is empty, read photons until the time gets out of the correlation window, and feed start buffer and the stop photons array in the process
        if (stop_buff_length == 0) {
            while (timetag < correlation_window_end && (*RecNum < RecNum_stop && *RecNum < NumRecords)) {
                next_photon(filehandle, record_type, RecNum, NumRecords, &oflcorrection, &timetag, &channel);
                // start photon -> store in start photon buffer (to be used later)
                if (channel == channel_start) {
                    push(start_buff_head, timetag, &start_buff_length);
                }
                // stop photon
                else {
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
        for(i = 0; i < stop_corr_buff_length; i++) {
            if(stop_corr_buff_head->next->val < start_time) {
                pop(stop_corr_buff_head, &stop_corr_buff_length);
            }
            else {
                break;
            }
        }
        // perform a histogram of the stop times - start time and add it to the main histogram result
        current = stop_corr_buff_head->next;
        while(current != NULL) {
            if (current->val - start_time < correlation_window) {
                i = (int) (current->val - start_time) * nb_of_bins / correlation_window;
                histogram[i] = histogram[i] + 1;
            }
            current = current->next;
        }
    }
}

void calculate_g2_fast(FILE* filehandle, long long record_type, int end_of_header, uint64_t *RecNum, uint64_t NumRecords, uint64_t RecNum_start, uint64_t RecNum_stop, uint64_t *time_vector, int *histogram, int nb_of_bins, int channel_start, int channel_stop)
{
    /*
     calculate_g2_fast() computes the g2 directly reading the measurement file. It uses a simple algorithm which stops at the first stop photon (histogram mode style). This algorithm is fast but loses some information and can exhibit an exponential decay artefact linked to the photon rates.
     Inputs:
     filehandle         FILE pointer with an open record file to read the photons
     record_type        record type which depends on the device which recorded the file (see constants at the beginning of file)
     end_of_header      offset in bytes to the beginning of the record section in the file
     RecNum             pointer to the index of the record being read
     NumRecords         total number of records
     RecNum_start       start of the section of records to analyse for the g2 (in terms of record index)
     RecNum_stop        stop of the section of records to analyse for the g2
     time_vector        precalculated array of times used for the x-axis of the g2 histogram. Should have nb_of_bins + 1 elements.
     histogram          preallocated array of zeros used for the g2 histogram. Should have nb_of_bins elements.
     nb_of_bins         number of bins for the histogram (should correspond to the length of the histogram array)
     channel_start      channel number used for start photons (sync will generally be 0)
     channel_stop       channel number used for stop photons (> 0, often 1)
     Outputs:
     filehandle         FILE pointer with reader at the position of last analysed record
     RecNum             index of last analysed record
     histogram          calculated g2 histogram
     */
    
    uint64_t start_time = 0;
    uint64_t stop_time = 0;
    uint64_t oflcorrection = 0;
    uint64_t timetag = 0;
    int channel = -1;
    int i = 0;
    uint64_t correlation_window = 0;
    //    long next_print = 0;
    correlation_window = time_vector[nb_of_bins];
    
    // reset file reader and go to the start position RecNum_start
    c_fseek(filehandle, end_of_header + 4 * RecNum_start);
    *RecNum = RecNum_start;
    
    // go to the start position RecNum_start
    while(*RecNum < RecNum_stop && *RecNum < NumRecords){
        //        if (*RecNum > next_print){
        //            printf("%ld/%ld\n", *RecNum, RecNum_stop);
        //            next_print = next_print + 1000000;
        //        }
        
        // FIND NEXT START PHOTON
        channel = -1;
        while((*RecNum < RecNum_stop && *RecNum < NumRecords) && channel != channel_start){
            next_photon(filehandle, record_type, RecNum, NumRecords, &oflcorrection, &timetag, &channel);
        }
        if (*RecNum >= RecNum_stop || *RecNum >= NumRecords){
            break;
        }
        // found a start photon
        else {
            start_time = timetag;
        }
        
        // FIND NEXT STOP PHOTON
        while ((*RecNum < RecNum_stop && *RecNum < NumRecords) && channel != channel_stop) {
            next_photon(filehandle, record_type, RecNum, NumRecords, &oflcorrection, &timetag, &channel);
        }
        // found a stop photon
        if (channel == channel_stop) {
            stop_time = timetag;
        }
        
        // ADD DELAY TO HISTOGRAM
        // add occurence to result histogram if the delay is in the correlation window
        if (stop_time - start_time < correlation_window) {
            i = (int) (stop_time - start_time) * nb_of_bins / correlation_window;
            histogram[i] = histogram[i] + 1;
        }
    }
}

