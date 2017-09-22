# some convenience tools
#

import logging
import numpy as np
import h5py
import platform
import os
import time
import datetime

try:
    import qt
    datadir = qt.config['datadir']
    #print datadir
except:
    # Added Mac compatibility. Does require data to be saved in correct folder (as below).
    # Added Linux compatibility, as well
    if os.name == 'posix':
        if (platform.system()=='Linux'):
            datadir = r'/home/'+os.environ.get('USER')+r'/Work/Data/'
        else:
            datadir = r'/Users/'+os.getlogin()+r'/Documents/teamdiamond/data'
    else:
        datadir = r'd:\measuring\data'
    #print datadir


    def get_file_list(some_path):
        orgin_path=os.getcwd()
        os.chdir(some_path)
        path = os.getcwd()
        file_list = [f for f in os.listdir('.') if (os.path.isfile(f) and f.endswith ('.pys'))]
        os.chdir(orgin_path)
        return file_list

    def find_string_in_array(f,seq):
      """Return items containing f in sequence """
      regex=re.compile(".*"+f+".*")#using regulary expressions for finding part of string in a list
      item=[m.group(0) for l in seq for m in [regex.search(l)] if m]
      return item[0] #item itself is a 1d-array so [0] has to be added

    def return_value(k, string):
        """Returns a wanted number in front of a string as an integer found in a given string k, output limited to 10 digits"""
        v_idx = k.find(string)
        #print 'string ',string,' found at position: ',v_idx
        returnvalue = 0
        for i in range(1,10):        
            try:
                returnvalue = float(k[v_idx-i:v_idx])
            except:
                pass
        return returnvalue

def nearest_idx(array, value):
    '''
    find the index of the value closest to the specified value.
    '''
    return np.abs(array-value).argmin()

def nearest_value(array, value):
    '''
    find the value in the array that is closest to the specified value.
    '''
    return array[nearest_idx(array,value)]

def get_timestamp_from_now():
    return timestamp_from_datetime(datetime.datetime.now())

def timestamp_from_datetime(datetime_):
    return datetime_.strftime('%Y%m%d%H%M%S')

def datetime_from_timestamp(timestamp):
    return datetime.datetime.strptime(timestamp,'%Y%m%d%H%M%S')
    
def verify_timestamp(timestamp):
    if len(timestamp) == 6:
        daystamp = time.strftime('%Y%m%d')
        tstamp = timestamp
    elif len(timestamp) == 14:
        daystamp = timestamp[:8]
        tstamp = timestamp[8:]
    elif len(timestamp) == 15: #### In case day and timestamp separted by _
        daystamp = timestamp[:8]
        tstamp = timestamp[9:]
    else:
        raise Exception("Cannot interpret timestamp '%s'" % timestamp)

    return daystamp, tstamp


def is_older(ts0, ts1):
    '''
    returns True if timestamp ts0 is an earlier data than timestamp ts1,
    False otherwise.
    '''
    if ts0 == None or ts1 == None:
        return True
    else:

        dstamp0, tstamp0 = verify_timestamp(ts0)
        dstamp1, tstamp1 = verify_timestamp(ts1)

        return (dstamp0+tstamp0) < (dstamp1+tstamp1)

def latest_data(contains='', older_than=None, newer_than=None,return_timestamp = False,raise_exc = True, folder=None, return_all=False, VERBOSE=False):
    '''
    finds the latest taken data with <contains> in its name.
    returns the full path of the data directory.

    if older_than is not None, than the latest data that fits and that
    is older than the date given by the timestamp older_than is returned.
    if newer_than is not None, than the latest data that fits and that
    is newer than the date given by the timestamp newer_than is returned

    If no fitting data is found, an exception is raised. Except when you specifically ask not to to
    this in: raise_exc = False, then a 'False' is returned.

    return_all = True: returns all the folders that satisfy the requirements (Cristian)
    '''

    if (folder==None):
        search_dir = datadir
    else:
        search_dir = folder

    daydirs = os.listdir(search_dir)

    if len(daydirs) == 0:
        logging.warning('No data found in datadir')
        return None

    daydirs.sort()

    # MAB 15-4-15 Added for weird file mac is automatically creating
    try:
        daydirs.remove('.DS_Store')
    except ValueError:
        pass

    measdirs = []
    valid_folders = []
    end_search = False
    i = len(daydirs)-1

    while i >= 0 and not end_search:  #Go through all days, starting from the latest one
        daydir = daydirs[i]
        if VERBOSE:
            print 'search dir: ', search_dir, ' day dir: ', daydir, ' break? ', end_search
        all_measdirs = [d for d in os.listdir(os.path.join(search_dir, daydir))]
        all_measdirs.sort() # get a sorted list of all measurements on that day

        for d in all_measdirs[::-1]:
            # this routine verifies that any output directory is a 'valid' directory
            # _timestamp = daydir + d[9:15] ### doing it like this breaks all other analysis scripts. Sorry cavities
            _timestamp = daydir + d[:6]
            try:
                dstamp,tstamp = verify_timestamp(_timestamp)
            except:
                if VERBOSE:
                    print 'Timestamp not valid: ', dstamp,tstamp
                continue
            timestamp = dstamp+tstamp

            if contains in d:    # file matching the name search found. Does it match the date constraints?
                if VERBOSE:
                    print 'contains ', contains, ' d ', d 
                    print ' break? ', end_search, ' return all? ', return_all           
                if older_than != None:
                    if not is_older(timestamp, older_than):
                        continue    # go to next measurement directory
                if newer_than != None:
                    if is_older(timestamp, newer_than):
                        if VERBOSE:
                            print 'Now I get to files which are too old'
                        end_search = True
                        continue
                measdirs.append(d) # if no continue command has been set, the file is valid.
                valid_folders.append(os.path.join(os.path.join(search_dir, daydir), d))
                if not return_all:
                    end_search = True
                if VERBOSE:
                    print 'Found a file. older than: ', older_than, ', timestamp: ', timestamp, ' Stop? ', end_search
            if end_search:
                break  # stop for loop
        i -= 1   #go to next day

    if VERBOSE:
        print 'measdirs length: ', len(measdirs)
    if len(measdirs) == 0:
        if raise_exc == True:
            raise Exception('No fitting data found containing {}.'.format(contains))
        else:
            return False
    else:
        measdirs.sort()
        if return_all:
            return valid_folders
            #return search_dir,daydir,measdirs
        measdir = measdirs[-1]
        if return_timestamp == False:
            return os.path.join(search_dir,daydir,measdir)
        else: 
            return str(daydir)+str(measdir[:6]) , os.path.join(search_dir,daydir,measdir)

def data_from_time(timestamp, folder = None):
    '''
    returns the full path of the data specified by its timestamp in the
    form YYYYmmddHHMMSS.
    '''
    global datadir

    if (folder != None):
        datadir = folder
        
    daydirs = os.listdir(datadir)

    if len(daydirs) == 0:
        raise Exception('No data in the data directory specified')

    daydirs.sort()
    daystamp, tstamp = verify_timestamp(timestamp)

    if not os.path.isdir(os.path.join(datadir,daystamp)):
        logging.warning("Requested day '%s' not found" % daystamp)
        return None

    measdirs = [ d for d in os.listdir(os.path.join(datadir,daystamp)) \
            if d[:6] == tstamp ]

    if len(measdirs) == 0:
        logging.warning("Requested data '%s'/'%s' not found" \
                % (daystamp, tstamp))
        return None
    elif len(measdirs) == 1:
        return os.path.join(datadir,daystamp,measdirs[0])
    else:
        logging.warning('Timestamp is not unique: ', str(measdirs))
        return None


def data_from_label(label, folder = None):
    
    if (folder != None):
        datadir = folder

    daydirs = os.listdir(datadir)

    if len(daydirs) == 0:
        raise Exception('No data in the data directory specified')

    daydirs.sort()
    
    if not os.path.isdir(os.path.join(datadir,daystamp)):
        logging.warning("Requested day '%s' not found" % daystamp)
        return None

    measdirs = [ d for d in os.listdir(os.path.join(datadir,daystamp)) \
            if d[:6] == tstamp ]

    if len(measdirs) == 0:
        logging.warning("Requested data '%s'/'%s' not found" \
                % (daystamp, tstamp))
        return None
    elif len(measdirs) == 1:
        return os.path.join(datadir,daystamp,measdirs[0])
    else:
        logging.warning('Timestamp is not unique: ', str(measdirs))
        return None

def measurement_filename(directory=os.getcwd(), ext='hdf5'):
    dirname = os.path.split(directory)[1]
    fn = dirname+'.'+ext

    if os.path.exists(os.path.join(directory,fn)):
        return os.path.join(directory,fn)
    else:
        logging.warning("Data path '%s' does not exist" % \
                os.path.join(directory,fn))
        return None

def get_date_time_string_from_folder(folder):
    if not os.path.isdir(folder):
        logging.error('Argument {} is not a folder'.format(folder))
    head,tail=os.path.split(folder)
    d = os.path.split(head)[1]
    if len(tail) < 6:
        logging.error('Argument {} is not a valid measurement folder'.format(folder))
    t = tail[:6]
    return d,t

def get_datetime_from_folder(folder):
    d,t = get_date_time_string_from_folder(folder)
    return datetime_from_timestamp(d+t)


def get_measurement_name_from_folder(folder):
    return os.path.split(folder)[1][7:]

def get_plot_title_from_folder(folder):
    d,t=get_date_time_string_from_folder(folder)
    timestamp = d + '/' + t
    measurementstring = get_measurement_name_from_folder(folder)
    default_plot_title = timestamp+'\n'+measurementstring
    return default_plot_title


def file_in_folder(folder, timestamp):
    all_files = [f for f in os.listdir(folder) if timestamp in f]    
    all_files.sort()

    if (len(all_files)>1):
        print 'More than one file satisfies the requirements! Import: '+all_files[0]

    return all_files[0]




############### 2014-06-11, Hannes: here i am putting some of wolgang's teleporation tools: file management etc. 

def get_all_msmt_filepaths(folder, suffix='hdf5', pattern=''):
    filepaths = []
    suffixlen = len(suffix)
    
    for root,dirs,files in os.walk(folder):
        for f in files:
            if len(f) > suffixlen and f[-suffixlen:] == suffix and pattern in f:
                filepaths.append(os.path.join(root, f))
                
    filepaths = sorted(filepaths)
    
    return filepaths

def get_msmt_name(pqf):
    """
    This assumes that there is only one group, whose name is the msmt name.
    """

    if type(pqf) == h5py._hl.files.File: 
        for k in pqf.keys():
            if pqf.get(k, getclass=True) == h5py._hl.group.Group and k in str(pqf):
                return k

        raise Exception('Cannot find the name of the measurement.')


    elif type(pqf) == str:
        _root, fn = os.path.split(pqf)

        f = h5py.File(pqf, 'r')

        for k in f.keys():
            if f.get(k, getclass=True) == h5py._hl.group.Group and k in fn:
                f.close()
                return k
    
        raise Exception('Cannot find the name of the measurement.')

    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise



def get_msmt_fp(folder, ext='hdf5'):
    dirname = os.path.split(folder)[1]
    fn = dirname+'.'+ext
    return os.path.join(folder, fn)

def get_msmt_header(fp):
    _root, fn = os.path.split(fp)
    root, folder = os.path.split(_root)
    daystamp = os.path.split(root)[1]
    return daystamp + '/' + folder

def has_analysis_data(pqf, name, analysisgrp = 'analysis', subgroup=None):
    """
    Returns a boolean stating if a file has analysis data or not. Enter a file or 
    a filepath and the name corresponding to the analysis data.
    """
    if type(pqf) == h5py._hl.files.File:
        if analysisgrp in pqf.keys():
            agrp = pqf.require_group(analysisgrp + ('/' + subgroup if subgroup!=None else ''))
            if name in agrp.keys():
                return True
            else:
                return False
        else:
            return False
    elif type(pqf) == str:
        try:
            f = h5py.File(pqf, 'r')
        except:
            print "Cannot open file", pqf
            raise
        
        if analysisgrp in f.keys():
            agrp = f.require_group(analysisgrp + ('/' + subgroup if subgroup!=None else ''))
            if name in agrp.keys():
                f.close()
                return True
            else:
                f.close()
                return False
        else:
            f.close()
            return False
    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise 

def get_analysis_data(pqf, name, analysisgrp = 'analysis', subgroup=None):
    """
    Return the analysis data from a file (1st entry filepath or file) withing 
    analysis with a specific name (2nd entry). Also returns the attributes of the file
    in an dictionary. 
    """
    if type(pqf) == h5py._hl.files.File:   
        agrp = pqf.require_group(analysisgrp + ('/' + subgroup if subgroup!=None else ''))

        if name not in agrp.keys():
            return None
        
        dat = agrp[name].value
        attrs = {}
        for (an, av) in agrp[name].attrs.items():
            attrs[an] = av
            print attrs
            print type(attrs)

        return dat, attrs
    elif type(pqf) == str:
        try:
            f = h5py.File(pqf, 'r')
        except:
            print "Cannot open file", pqf
            raise
        
        agrp = f.require_group(analysisgrp + ('/' + subgroup if subgroup!=None else ''))

        if name not in agrp.keys():
            return None
    
        dat = agrp[name].value
        attrs = {}
        for (an, av) in agrp[name].attrs.items():
            attrs[an] = av
        
        f.close()
        
        return dat, attrs
    else:
        print "Neither filepath nor file enetered in function please check:", pqf
        raise


def clear_analysis_data(fp, ANALYSISGRP = 'analysis'):
    """
    Deletes analysis group from the hdf5 file. Can be used 
    to delete any group by changing the name put in the function.
    """
    try:
        f = h5py.File(fp, 'r+')
    except:
        print "Cannot open file", fp
        raise

    if ANALYSISGRP in f.keys():
        try:
            del f[ANALYSISGRP]
            f.flush()
            f.close()
        except:
            f.close()
            raise
    else:
        f.close()

def clear_raw_data(fp, name):
    """
    Deletes the raw data in the main folder. 
    """

    try:
        f = h5py.File(fp, 'r+')
    except:
        print "Cannot open file", fp
        raise

    if name in f.keys():
        try:
            del f[name]
            f.flush()
            f.close()
        except:
            f.close()
            raise
    else:
        f.close()

def set_raw_data(fp, name, data):
    """
    Save the data in the main folder. Does not save any attributes
    """
    try:
        f = h5py.File(fp, 'r+')
    except:
        print "Cannot open file", fp
        raise

    f[name] = data
    f.flush()
    f.close()        
    
def set_analysis_data(fp, name, data, attributes, subgroup=None, ANALYSISGRP = 'analysis', permissions='r+'):
    """
    Save the data in a subgroup which is set to analysis by default and caries the name
    put in the function. Also saves the attributes.
    """
    try:
        f = h5py.File(fp, permissions)
    except:
        print "Cannot open file", fp
        raise

    try:
        agrp = f.require_group(ANALYSISGRP + ('/' + subgroup if subgroup!=None else ''))

        
        if name in agrp.keys():
            del agrp[name]
        agrp[name] = data
        f.flush()
        
        for k in attributes:
            agrp[name].attrs[k] = attributes[k]
            #print agrp[name].attrs[k]
        #    agrp[name].attrs[k] = kw[k]
        
        f.flush()
        f.close()        
    except:
        f.close()
        raise

def validate_ip(s):
    a = s.split('.')
    if len(a) != 4:
        return False
    for x in a:
        if not x.isdigit():
            return False
        i = int(x)
        if i < 0 or i > 255:
            return False
    return True