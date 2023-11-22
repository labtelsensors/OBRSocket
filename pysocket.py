import socket 
import time
import matplotlib.pyplot as plt
from pathlib import Path
import logging
import pandas as pd
import numpy as np
import os
import math

# Server parameters
HOST = '192.168.0.199'
PORT = 1

# Set a timeout for blocking socket operations
socket.timeout(10)

# Msg maximum size
msg_size = 1024 * 50 * 2

# Code parameters
file_dir = "C:/Users/Mariana Lyra/Documents/Luna/data/cao/strain/distributed_strain5"
file_name = "0d.txt"
move_trace = False
measure_scan = True
save_figure = True
fig_names = ["freq_domain.png","time_domain_amplitude.png","time_domain_spectralshift.png"]

# OBR Dictionary
domain_dict = {0:"Time", 1: "Frequency"}

type_dict = {   'Time':         {'x':  {0: "Time (ns)", 1: "Length (m)", 2: "Length (ft)",
                                        3: "Length (in),", 4: "Length (mm)"
                                        },
                                 'y':   {0: "Amplitude (dB/mm)", 1: "Linear Amplitude", 2: "Polarization States",
                                        3: "Phase Derivative (nm)", 4: "Amplitude (dB)", 5: "Magnitude Difference (db/mm)",
                                        6: "Spectral Shift (GHz)", 7: "Spectral Shift Quality"
                                        }
                                },
                'Frequency':    {'x':   {0: "Wavelength (nm)", 1: "Frequency (GHz)", 2: "Frequency (THz)"
                                        },
                                'y':    {0: "Return Loss (dB)", 1: "Linear Amplitude (1/mm)", 2: "Polarization States (dB)",
                                        3: "Group Delay (ns)"
                                        } 
                                }

            }

# Maximum trials to complete scan
max_trial = 200
trial = 0
opc = 0

# Trace sweeping parameters
scan_start = 2.75
scan_end = 3.25
scan_div = 0.5

# Minimum frequency amplitude
min_freq_amplitude = -100.00

def getQuotientRemainder(n, divisor):
    res =  n/divisor
    quo = round(res)
    rem = 0 if abs(quo-res) < np.finfo(float).resolution else res-quo
    # Update quotient because the total number of scans is (quotient + 1)
    # quo += 1
    return (quo,rem)

def sendCommand(s,msg,msg_size=msg_size):
    '''
    Input:
        - s: socket - socket instance
        - msg: str - command
        - msg_size - integer - mesage size
    Returns:
        - data - str - received msg
    '''
    s.sendall(msg)
    time.sleep(0.2)
    if (b"?") in msg:
        data = s.recv(msg_size)
        return data
    else:
        return

def byte2floatList(msg):
    data = list()
    # Split data for each line terminator
    _msg = msg.split(b'\r')  
    
    # Split data for each line separator
    _msg = [line.split(b'\t') for line in _msg]

    # Convert to string
    _msg = [[str(item,'utf-8') for item in line] for line in _msg]

    # Remove null string terminator
    _msg = [[item.replace('\x00','') if('\x00' in item) else item for item in line] for line in _msg]    
    
    for line in _msg:
        try:
            line_data = [float(item) for item in line]
            data.append(line_data)
        except:
            for i, item in enumerate(line):
                # Look for msg terminator
                if ('///' in item):
                    logging.info("End of message detected: " + item)
                # Removes message that appers in frequency domain, when data is below the minimum OBR limit
                elif ('#IND0000' in item): 
                    data.append([min_freq_amplitude])
                    logging.info("Could not parse msg '{}' to data. Appending the minimum frequency amplitude instead".format(item))
                else:    
                    # If end of line
                    if (i == len(line)-1):
                        logging.info("Could not parse msg '{}' to data".format(line))
    return data

def list2dataFrame(x, y, labels):
    x = np.array(x)
    x_y = x
    x_y = np.concatenate((x_y,np.array(y)),axis=1)
    df = pd.DataFrame(x_y, columns=labels)
    return df

def getInfoFromCommand(command):
    # Check command
    axis = ''
    if 'FETC:XAXI?' in command:
        axis = 'x'
    elif 'FETC:MEAS?' in command:
        axis = 'y'
    elif 'FETC:RAW?' in command:
        return dict([('trace', command[10]), ('label', 'Polarization States')])
    else:
        logging.error("Invalid command.")
        return
    
    # Extract information
    try:
        trace, domain, graph_type, dec = command.split('?')[-1].split(',')
    except ValueError:
        logging.error("Information could not be extracted from command.")
        return
    
    # Convert to label
    domain = domain_dict[int(domain)]
    graph_type = type_dict[domain][axis][int(graph_type)]  

    # Return dictionary
    info = [('trace', trace), ('domain', domain), ('label', graph_type), ('dec', dec)]
    return dict(info)   

def createDir(dir, file_name, folder=[]):
    if folder:
        dir = os.path.join(dir,folder)
    path = os.path.join(dir,file_name)
    return path   

def saveFile(path, df, header = [], float_format = "%.8f"):
    '''
    Input:
        path - string - path for saving file
        df - dataframe - data to be saved
    '''
    filepath = Path(path)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(path,'w', newline='') as file:
        if header:
            file.write(header + '\n')
        df.to_csv(file,sep='\t',index=False,float_format=float_format)
    return

def plotFig(x, y, labels, fig_path):
    if (len(x)!=len(y)) & (len(x)!=len(labels)) & (len(x)!=len(fig_path)):
        raise Exception("Input parameters must have the same sime.")
    
    for i in range(len(x)):
        figpath = Path(fig_path[i])
        figpath.parent.mkdir(parents=True, exist_ok=True)

        fig = plt.figure(1)
        plt.plot(x[i], y[i],'-') 
        plt.xlabel(labels[i][0])
        plt.ylabel(labels[i][1])
        plt.grid(alpha=0.5,linestyle='--')
        plt.savefig(figpath)
        #plt.show()
        plt.close(fig)
    return

def xArraySpectralShift(start, range, ndiv):
    """ Interpolate data to get measurement in a specific length.
    Args:
        df - dataframe
    Returns:
        ynew - numpy array - interpolated measurement
    """
    x = np.linspace(start, start+range, ndiv)
    x = [[item] for item in x]
    return x


# Check if parameters are consistent
if (scan_start > scan_end):
    raise Exception("The sweeping end ({} m) must be greater than the sweeping start ({}} m).".format(scan_end, scan_start))

# Calculate remaining scanning parameters
n_scans, n_rem = getQuotientRemainder((scan_end-scan_start), scan_div)
print("n_scans: {}| n_rem: {}".format(n_scans,n_rem))

if (n_rem != 0):
    logging.info("The sweeping division length is not divible by the sweeping interval. The sweeping will end at {}".format(scan_start+scan_div*n_scans))

# Connects to server
with socket.socket(socket.AF_INET,socket.SOCK_STREAM) as s:
    s.connect((HOST,PORT))
    
    # Send command: System Version
    commnd = "SYST:VER?"
    data = sendCommand(s, commnd.encode())
    print("{}".format(data.decode()))

    # Enable Sensing # CONF:SENS 1
    _ = sendCommand(s, b"CONF:SENS 1")

    # Sensing Range Unit # CONF:XUNI
    _ = sendCommand(s, b"CONF:XUNI 1")

    # Sensing Range # CONF:INTW
    commnd = "CONF:INTW {:.3f}".format(scan_div)
    _ = sendCommand(s, commnd.encode())

    # Genge length CONF:GLEN
    _ = sendCommand(s, b"CONF:GLEN 2.0")

    # Sensor spacing CONF:SSPA
    _ = sendCommand(s, b"CONF:SSPA 2.0")
    
    if measure_scan:
        # Perform SCAN
        time.sleep(0.2)
        commnd = "SCAN"
        _ = sendCommand(s, commnd.encode())
        # Wait until device is ready
        while((opc==0) & (trial < max_trial)):
            trial += 1
            opc = int(sendCommand(s, "*OPC?".encode()).decode().replace('\x00',''))

    # Move trace from A to D
    if move_trace:
        commnd = "SYST:COPY A, D"
        _ = sendCommand(s, commnd.encode())
        print('Moving to Reference Trace.')

    # Data Processing
    for i in range(n_scans):        
        # Cursor centralization CONF:INTC
        cursor_start = scan_start + scan_div * i
        print("Cursor starts at {:.3f} m".format(cursor_start))
        commnd = "CONF:INTS {:.3f}".format(cursor_start)
        _ = sendCommand(s, commnd.encode())
 
        # Send command: Get X measuremen for frequency domain
        commnd = "FETC:XAXI? A, 1, 0, 1"
        x_freq = sendCommand(s, commnd.encode())
        x_freq = byte2floatList(x_freq)
        info_x =  getInfoFromCommand(commnd)

        # Send command: Get Y measurement for frequency domain
        commnd = "FETC:MEAS? A, 1, 0, 1"
        y_freq = sendCommand(s, commnd.encode(),msg_size)
        y_freq = byte2floatList(y_freq)
        info_y =  getInfoFromCommand(commnd)

        labels_freq = [info_x['label'], info_y['label']]
        df_freq = list2dataFrame(x_freq, y_freq, labels_freq) 

        # Send command: Get X measuremen for time domain
        commnd = "FETC:XAXI? A, 0, 1, 1"
        x_time = sendCommand(s, commnd.encode())
        x_time = byte2floatList(x_time)
        info_x =  getInfoFromCommand(commnd)

        # Send command: Get Y measurement for time domain
        commnd = "FETC:MEAS? A, 0, 0, 1"
        y_time = sendCommand(s, commnd.encode(),msg_size)
        y_time = byte2floatList(y_time)
        info_y =  getInfoFromCommand(commnd)

        labels_time = [info_x['label'], info_y['label']]
        df_time = list2dataFrame(x_time, y_time, labels_time) 

        # Send command: Get Y measurement for time domain - spectral shift
        commnd = "FETC:MEAS? A, 0, 6, 1"
        y_time_spec = sendCommand(s, commnd.encode(),msg_size)
        y_time_spec = byte2floatList(y_time_spec)
        info_y =  getInfoFromCommand(commnd)

        # Send command: Get X measuremen for time domain - spectral shift
        x_time_spec = xArraySpectralShift(cursor_start,scan_div,len(y_time_spec))
        labels_time_spec = [type_dict['Time']['x'][1], info_y['label']]
        df_time_spec = list2dataFrame(x_time_spec, y_time_spec, labels_time_spec) 

        # Concatenate dataframes and lists
        df = pd.concat([df_freq, df_time, df_time_spec],axis=1)
        #df = pd.concat([df, df_time_spec], axis=1)
        labels = [labels_freq, labels_time, labels_time_spec]
        x = [x_freq, x_time, x_time_spec]
        y = [y_freq, y_time, y_time_spec]

        # Get file hearder
        header = sendCommand(s, "FETC:MDET?".encode()).decode()

        # Save file
        file_path = createDir(file_dir, file_name, folder="{:.3f}m".format(cursor_start))
        saveFile(file_path, df, header)
        print("File successfully saved in '{}'.".format(file_path))
        
        if save_figure:
            fig_path = [createDir(file_dir, fig_name, folder="{:.3f}m".format(cursor_start))
                        for fig_name in fig_names]
            plotFig(x, y, labels, fig_path)

    # Send command: Quit
    s.sendall(b"*quit")
    time.sleep(0.2)
    s.close()

