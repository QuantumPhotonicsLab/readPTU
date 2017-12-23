#ifndef PARSERS_C_
#define PARSERS_C_

static inline void ParsePHT2(uint32_t record, int * channel, uint64_t * timetag, uint64_t * oflcorrection);
static inline void ParseHHT2_HH1(uint32_t record, int * channel, uint64_t * timetag, uint64_t * oflcorrection);
static inline void ParseHHT2_HH2(uint32_t record, int * channel, uint64_t * timetag, uint64_t * oflcorrection);

static inline void ParsePHT2(uint32_t record, int * channel,
                             uint64_t * timetag, uint64_t * oflcorrection)
{
    /*
     ProcessPHT2() reads the next records of a file until it finds a photon, and then returns.
     Inputs:
     filehandle         FILE pointer with an open record file to read the photons
     oflcorrection      pointer to an unsigned integer 64 bits. Will record the time correction in the timetags due to overflow (see output for details).
     buffer             buffer from which to read the next record (file chunk read buffer)
     Outputs:
     filehandle         FILE pointer with reader at the position of last analysed record
     oflcorrection      offset time on the timetags read in the file, due to overflows.
     buffer             buffer of the next chunk of records, containing for each a timetag and a channel number.
     If a photon is read, timetag of this photon. Otherwise, timetag == 0. It already includes the overflow correction
     so the value can be used directly.
     If a photon is read, channel of this photon. 0 will usually be sync and >= 1 other input channels. If the record is
     not a photon, channel == -1 for an overflow record, -2 for a marker record.
     */
    /* FUNCTION TESTED QUICKLY */
    
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
    
    Record.allbits = record;
    
    if(Record.bits.channel == 0xF) //this means we have a special record
    {
        //in a special record the lower 4 bits of time are marker bits
        markers = Record.bits.time & 0xF;
        if(markers == 0) //this means we have an overflow record
        {
            *timetag = 0;
            *channel = -1;
            *oflcorrection += T2WRAPAROUND; // unwrap the time tag overflow
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

static inline void ParseHHT2_HH1(uint32_t record, int * channel,
                                 uint64_t * timetag, uint64_t * oflcorrection)
{
    /*
     ProcessHHT2() reads the next records of a file until it finds a photon, and then returns.
     Inputs:
     filehandle         FILE pointer with an open record file to read the photons
     HHVersion          Hydrahard version 1 or 2. Depends on record type specification in the header.
     oflcorrection      pointer to an unsigned integer 64 bits. Will record the time correction in the timetags due to overflow (see output for details).
     buffer             buffer from which to read the next record (file chunk read buffer)
     Outputs:
     filehandle         FILE pointer with reader at the position of last analysed record
     oflcorrection      offset time on the timetags read in the file, due to overflows.
     buffer             buffer of the next chunk of records, containing for each a timetag and a channel number.
     If a photon is read, timetag of this photon. Otherwise, timetag == 0. It already includes the overflow correction
     so the value can be used directly.
     If a photon is read, channel of this photon. 0 will usually be sync and >= 1 other input channels. If the record is
     not a photon, channel == -1 for an overflow record, -2 for a marker record.
     */
     //FUNCTION TESTED 
    
    const uint64_t T2WRAPAROUND_V1 = 33552000;
    union{
        uint32_t   allbits;
        struct{ unsigned timetag  :25;
            unsigned channel  :6;
            unsigned special  :1; // or sync, if channel==0
        } bits;
    } T2Rec;
    
    T2Rec.allbits = record;
    
    unsigned ch = T2Rec.bits.channel;
    unsigned sp = T2Rec.bits.special;
    unsigned tm = T2Rec.bits.timetag;

    // Negative Channels represent overflows or markers
    *oflcorrection += T2WRAPAROUND_V1 * (ch==0x3F);
    *channel = (!sp) * (ch + 1) - sp * ch;
    *timetag = *oflcorrection + tm;
}

static inline void ParseHHT2_HH2(uint32_t record, int *channel,
                                 uint64_t *timetag, uint64_t *oflcorrection)
{
    /*
     ProcessHHT2() reads the next records of a file until it finds a photon, and then returns.
     Inputs:
         filehandle         FILE pointer with an open record file to read the photons
         HHVersion          Hydrahard version 1 or 2. Depends on record type specification in the header.
         oflcorrection      pointer to an unsigned integer 64 bits. Will record the time correction in the timetags due to overflow (see output for details).
         buffer             buffer from which to read the next record (file chunk read buffer)
     Outputs:
         filehandle         FILE pointer with reader at the position of last analysed record
         oflcorrection      offset time on the timetags read in the file, due to overflows.
         buffer             buffer of the next chunk of records, containing for each a timetag and a channel number.
                            If a photon is read, timetag of this photon. Otherwise, timetag == 0. It already includes the overflow correction 
                                so the value can be used directly.
                            If a photon is read, channel of this photon. 0 will usually be sync and >= 1 other input channels. If the record is 
                                not a photon, channel == -1 for an overflow record, -2 for a marker record.
     */
    /* FUNCTION TESTED */

    const uint64_t T2WRAPAROUND_V2 = 33554432;
    union{
        uint32_t   allbits;
        struct{ unsigned timetag  :25;
            unsigned channel  :6;
            unsigned special  :1;
        } bits;
    } T2Rec;

    T2Rec.allbits = record;
    
    unsigned ch = T2Rec.bits.channel;
    unsigned sp = T2Rec.bits.special;
    unsigned tm = T2Rec.bits.timetag;

    // Negative Channels represent overflows or markers
    // if(ch==0x3F) {  //an overflow record
        // if(T2Rec.bits.timetag!=0) {
                // *oflcorrection += T2WRAPAROUND_V2 * T2Rec.bits.timetag * (ch==0x3F);
            // }
            // else {  // if it is zero it is an old style single overflow
                // *oflcorrection += T2WRAPAROUND_V2;  //should never happen with new Firmware!
            // }
    // }
    *oflcorrection += T2WRAPAROUND_V2 * tm * (ch==0x3F);
    *channel = ch + 1 - sp*(2*ch + 1);
    *timetag = *oflcorrection + tm;
}

#endif /* PARSERS_C_ */