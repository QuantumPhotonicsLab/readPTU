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

// How big the file chunking will be
#define RECORD_CHUNK = 512;

// ================================================
// Buffer for keeping track of records
// ================================================
typedef struct {
    uint64_t timetag;
    int channel;
} photon;

typedef struct {
    // photon photons[RECORD_CHUNK]; // using this will improve memory locality
    uint64_t timetag[RECORD_CHUNK];
    int channel[RECORD_CHUNK];
    size_t head;  // to keep track of what was the last read record in buffer
    size_t count; // if don't have enough photons, for example due to many records being oflcorrection flags
} photon_buf_t;

void photon_buf_reset(photon_buf_t *buffer) {
    buffer->head = 0;
    buffer->count = 0;
}

void photon_buf_pop(photon_buf_t *buffer, uint64_t *timetag, int *channel) {
    size_t head = buffer->head;
//    *oflcorrection = buffer->oflcorrection[head];
    *timetag = buffer->timetag[head];
    *channel = buffer->channel[head];
    buffer->head = head + 1;
}

void photon_buf_push(photon_buf_t *buffer, uint64_t timetag, int channel) {
    size_t count = buffer->count;
    buffer->timetag[count] = timetag;
    buffer->channel[count] = channel;
    buffer->count = count + 1;
}
// ================================================
// END Buffer for keeping track of records
// ================================================

int c_fseek(FILE *filehandle, long int offset)
{
    return fseek(filehandle, offset, SEEK_SET);
}

void ProcessPHT2(FILE* filehandle, photon_buf_t *buffer,  uint64_t *oflcorrection,)
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
    uint32_t TTTRRecord[RECORD_CHUNK];
    
    fread(&TTTRRecord, RECORD_CHUNK, sizeof(TTTRRecord) ,filehandle);
    for(size_t i = 0; i < RECORD_CHUNK; i++) {
        Record.allbits = TTTRRecord[i];
    
        if(Record.bits.channel == 0xF) //this means we have a special record
        {
            //in a special record the lower 4 bits of time are marker bits
            markers = Record.bits.time & 0xF;
            if(markers == 0) //this means we have an overflow record
            {
                *oflcorrection += T2WRAPAROUND; // unwrap the time tag overflow
            }
        }
        else
        {
            if((int)Record.bits.channel > 4) //Should not occur
            {
//                buffer->timetag[i] = 0;                                                   ///
//                buffer->channel[i] = -3;                                                  ///
            }
            else
            {
                photon_buf_push(buffer,
                                *oflcorrection + Record.bits.time,
                                Record.bits.channel);
            }
        }
    }
}

void ProcessHHT2(FILE* filehandle, int HHVersion, photon_buf_t *buffer,  uint64_t *oflcorrection,)
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
    uint32_t TTTRRecord[RECORD_CHUNK];
    
    
    fread(&TTTRRecord, RECORD_CHUNK, sizeof(TTTRRecord) ,filehandle);
    
    for(size_t i = 0; i < RECORD_CHUNK; i++) {
        T2Rec.allbits = TTTRRecord[i];
    
        if(T2Rec.bits.special==1)
        {
            if(T2Rec.bits.channel==0x3F) //an overflow record
            {
                if(HHVersion == 1)
                {
                    *oflcorrection += T2WRAPAROUND_V1;
                }
                else
                {
                    //number of overflows is stored in timetag
                    if(T2Rec.bits.timetag==0) //if it is zero it is an old style single overflow
                    {
                        *oflcorrection += T2WRAPAROUND_V2;  //should never happen with new Firmware! ///
                    }
                    else
                    {
                        *oflcorrection += T2WRAPAROUND_V2 * T2Rec.bits.timetag; ///
                    }
                }
            }
            
            if((T2Rec.bits.channel>=1)&&(T2Rec.bits.channel<=15)) //markers
            {
//                buffer->timetag[i] = *oflcorrection + T2Rec.bits.timetag;  ///
                //Note that actual marker tagging accuracy is only some ns.
//                buffer->channel[i] = -2;   ///
                // *channel = T2Rec.bits.channel;
            }
            
            else if(T2Rec.bits.channel==0) //sync
            {
//                buffer->timetag[i] = *oflcorrection + T2Rec.bits.timetag;  ///
//                buffer->channel[i] = T2Rec.bits.channel;        ///
            }
        }
        else //regular input channel
        {
            photon_buf_push(buffer,
                            *oflcorrection + Record.bits.time,
                            Record.bits.channel);
    }
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


int next_photon(FILE* filehandle, long long record_type, uint64_t *RecNum, uint64_t NumRecords, photon_buf_t *buffer, uint64_t *oflcorrection, uint64_t *timetag, int *channel)
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
    
    // We will sacrifice up to 512 records at the end of the file in order to simplify the logic of the function.
    
    if (buffer->head < buffer->count) {
        // pop a photon
        photon_buf_pop(buffer, timetag, channel);
        return 1;
    } else {
        // we need to replenish the photon pool
        photon_buf_reset(buffer);
        while (buffer->count <= 0 && (*RecNum+RECORD_CHUNK) < NumRecords) {
            *RecNum = *RecNum + RECORD_CHUNK;
            
            switch (record_type) {
                case rtPicoHarpT2:
                    ProcessPHT2(filehandle, buffer);
                    break;
                case rtPicoHarpT3:
                    //ProcessPHT3(TTTRRecord);
                    break;
                case rtHydraHarpT2:
                    ProcessHHT2(filehandle, 1, buffer, oflcorrection);
                    *RecNum = *RecNum + 1;
                    break;
                case rtHydraHarpT3:
                    //ProcessHHT3(TTTRRecord, 1);
                    break;
                case rtHydraHarp2T2:
                case rtTimeHarp260NT2:
                case rtTimeHarp260PT2:
                    ProcessHHT2(filehandle, 2, buffer, oflcorrection);
                    break;
                case rtHydraHarp2T3:
                case rtTimeHarp260NT3:
                case rtTimeHarp260PT3:
                    //ProcessHHT3(TTTRRecord, 2);
                    break;
                default:
                    return 0;
            }
            
        };
        
        // Check if the last while loop added some photons and pop if possible
        if (buffer->count > 0) {
            photon_buf_pop(buffer, timetag, channel);
            return 1;
        } else {
            return 0;
        }
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
    photon_buf_t photon_buffer;
    photon_buf_reset(&photon_buffer);
    
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
            while(*RecNum < NumRecords) // this has to be change now we are looping as long as we get photons from next_photon
            {
                next_photon(filehandle, record_type, RecNum, NumRecords, &photon_buffer, &oflcorrection, &timetag, &channel));
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
