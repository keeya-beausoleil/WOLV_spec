
#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal # a signal processing module in the scipy package
import obspy
from obspy import read, read_inventory
import os
#%%
# data_path = 'C:\\Users\\dbehr\\Documents\\Seismograph_scripts\\Test_files\\'
# filename0 = 'AM.R7F0A.00.EHZ.D.2019.280'
# filename1 = 'AM.R7F0A.00.EHZ.D.2019.281'
# filename2 = 'AM.R7F0A.00.EHZ.D.2019.282'
# filename3 = 'AM.R7F0A.00.EHZ.D.2019.283'

data_path = '/data/stor/basic_data/seismic_data/day_vols/WOLVERINE/WOLN/'
file_prefix = 'WOLN.XX.00.HHZ.2022.' 
days = np.arange(120,125,1) 
stn_names = []
for day in days: 
    stn = data_path + file_prefix + str(day)
    if os.path.exists(stn):
        stn_names.append(stn)
    else:
        continue



'''
# filename0 = 'BBGU.LM..HHZ.2017.187'
#filename1 = 'SELC.XX..HHZ.2018.206'
#filename2 = 'SELC.XX..HHZ.2018.207'

data_path = '/Users/timb/Documents/syncs/OneDrive - University of Idaho/RESEARCHs/MoVE_Gulley_Greenland/Seismic_Data/SE53_early_june/'
filename0 = 'SE53.XX..HHZ.2018.156'

data_path = '/Users/timb/Downloads/'
filename0 = 'SE7.YG.00.HHZ.2021.225'
filenames = ['SE7.YG.00.HHZ.2021.223',
             'SE7.YG.00.HHZ.2021.224',
             'SE7.YG.00.HHZ.2021.225',
             'SE7.YG.00.HHZ.2021.226',
             'SE7.YG.00.HHZ.2021.227',
]
'''
#%%

st = obspy.Stream()
for stnname in stn_names:
    st += obspy.read(stnname)
# st += obspy.read(data_path + filename1)
# st += obspy.read(data_path + filename2)
# st += obspy.read(data_path + filename3)
# st += obspy.read('C:\\Users\\dbehr\\Documents\\Seismograph_scripts\\Test_files\\AM.R7F0A.00.EHZ.D.2019.284')
# st += obspy.read('C:\\Users\\dbehr\\Documents\\Seismograph_scripts\\Test_files\\AM.R7F0A.00.EHZ.D.2019.285')
# st += obspy.read('C:\\Users\\dbehr\\Documents\\Seismograph_scripts\\Test_files\\AM.R7F0A.00.EHZ.D.2019.286')
# st += obspy.read('C:\\Users\\dbehr\\Documents\\Seismograph_scripts\\Test_files\\AM.R7F0A.00.EHZ.D.2019.287')
# st += obspy.read('C:\\Users\\dbehr\\Documents\\Seismograph_scripts\\Test_files\\AM.R7F0A.00.EHZ.D.2019.288')

#%%
# sort
# st.sort(['starttime'])
#st.merge()

# start time in plot equals 0
# dt = st[0].stats.starttime.timestamp
#%%
# # Go through the stream object, determine time range in julian seconds
# # and plot the data with a shared x axis
# plt.figure()
# ax = plt.subplot(4, 1, 1)  # dummy for tying axis
# for i in range(3):
#     plt.subplot(4, 1, i + 1, sharex=ax)
#     t = np.linspace(st[i].stats.starttime.timestamp - dt,
#                     st[i].stats.endtime.timestamp - dt,
#                     st[i].stats.npts)
#     plt.plot(t, st[i].data)

# # Merge the data together and show plot in a similar way
# st.merge(method=1)
# plt.subplot(4, 1, 4, sharex=ax)
# t = np.linspace(st[0].stats.starttime.timestamp - dt,
#                 st[0].stats.endtime.timestamp - dt,
#                 st[0].stats.npts)
# plt.plot(t, st[0].data, 'r')
# plt.show()

#%%
pre_filt1 = [0.05, 0.1, 80, 100] #standard bandpass parameters to remove low and high freq. signal 
instr_resp = read_inventory("/data/stor/basic_data/seismic_data/day_vols/WOLVERINE/resp/Wolverine_station_20240824.xml")
station_rem=st.copy()
station_rem.remove_response(inventory=instr_resp, water_level= None,pre_filt=pre_filt1)

#%%
start = station_rem[0].stats.starttime
end = station_rem[0].stats.endtime
station_rem_spec = station_rem.copy()
#station_rem_cut.trim(starttime=start+(20*60*60), endtime = start+(20*60*60)+60)
#station_rem_spec = station_rem_cut.copy()
#print(station_rem_cut)
#%% Calculate the spectra of a segment of data.
num_welch_segs = 8

freqs, power = signal.welch(station_rem_spec[0].data, fs=station_rem_spec[0].stats.sampling_rate, window='hann', 
             nperseg=station_rem_spec[0].stats.npts/num_welch_segs, noverlap=None, nfft=None, 
             detrend='linear', return_onesided=True, 
             scaling='density', axis=-1, average='mean')

fig, ax = plt.subplots(num=3, clear=True)
station_rem_spec.plot(fig=fig)

ax_spec = fig.axes[0]
ax_spec.set_position([.1, .1, .85, .6 ])
wv_spec = fig.axes[1]
wv_spec.set_position([.1, .8, .85, .15 ])
#freqs, power = signal.periodogram(station_rem_spec[0].data, fs=station_rem_spec[0].stats.sampling_rate, window='hanning')
power_db = 10 * np.log10(power) # dB
ax_spec.plot(freqs, power_db,'r')
ax_spec.set_xscale('log')
ax_spec.set_xlabel('Frequency [Hz]')
ax_spec.set_ylabel('Power [dB ref. counts^2/Hz]')
#ax_spec.set_ylim(-40, 60)
ax_spec.grid()


#%%
'''
tr = st[0].copy()
# tr.plot()
starttime=obspy.UTCDateTime('2017-07-07T04:56:13.0')
starttime=obspy.UTCDateTime('2017-07-07T04:57:08.2')
starttime=obspy.UTCDateTime('2017-07-07T04:57:29.6') #surface
starttime=obspy.UTCDateTime('2017-07-07T04:58:01.3')
starttime=obspy.UTCDateTime('2017-07-07T05:08:11.7')
timespan = 1.0 # seconds


starttime=obspy.UTCDateTime('2017-07-06T06:20')
starttime=obspy.UTCDateTime('2017-07-07T01:00:00')

timespan = 10*60*60.0 # seconds

starttime=obspy.UTCDateTime('2018-06-05T12:00:00')


starttime=obspy.UTCDateTime('2021-08-13T02:00:00')


# starttime=obspy.UTCDateTime('2018-07-25T20:10')

# Trim or no?
# tr.trim(starttime=starttime, endtime=starttime + timespan) # Trim data

#%% Spectrogram with waveform too.
fig, ax = plt.subplots(num=3, clear=True, figsize=(16,10))
tr.plot(fig=fig) # Plot waveform into the figure fig

# Define and reshape axes
ax_spec = fig.axes[0] 
ax_spec.set_position([.1, .1, .85, .6 ])

wv_spec = fig.axes[1]
wv_spec.set_position([.1, .8, .85, .15 ])
wv_spec.set_ylim(-15000, 15000)
# wv_spec.set_ylim(-4000, 4000) # change these waveform axes limits to suit!

# freqs, power = signal.periodogram(data, fs=sample_rate, window='hanning', 
                                  # nfft=NFFT)
# NFFT is the number of data points to use in calculating periodograms.
#   It's usually a power of 2.  Smaller numbers provide better temporal resolution
#   and larger numbers provide better (low frequency) spectral resolution.
NFFT = 2**5
NFFT = 2**14
# Calculate the spectrogram.  vmin and vmax define the power bounds of the colormap
ax_spec.specgram(tr.data, NFFT=NFFT, Fs=int(tr.stats.sampling_rate), Fc=None, 
                 detrend='linear', window=None, noverlap=NFFT*.5, cmap=None, 
                 xextent=None, pad_to=None, sides=None, scale_by_freq=True, 
                 # mode='psd', scale=None, vmin=None, vmax=-40)
                #  mode='psd', scale=None, vmin=-15, vmax=45)
                 mode='psd', scale=None, vmin=-20, vmax=30)
ax_spec.set_ylabel('Frequency [Hz]')
ax_spec.set_xlabel('Time [s]')
# ax_spec.set_ylabel('Power [dB ref. counts^2/Hz]')
# ax_spec.set_xlim(0, timespan)

# ax_spec.set_yscale('log')
# ax_spec.set_ylim(0.5, 250)
# ax_spec.set_ylim(0.05, 125)
ax_spec.set_ylim(0, 25)
ax_spec.grid()

# ax_spec.set_ylim(0, 4)
# ax_spec.grid()

#%% Calculate the spectra of a segment of data.
num_welch_segs = 8

freqs, power = signal.welch(tr.data, fs=tr.stats.sampling_rate, window='hann', 
             nperseg=tr.stats.npts/num_welch_segs, noverlap=None, nfft=None, 
             detrend='linear', return_onesided=True, 
             scaling='density', axis=-1, average='mean')
power_db = 10 * np.log10(power) # dB


fig, ax = plt.subplots(num=3, clear=True)
tr.plot(fig=fig)

ax_spec = fig.axes[0]
ax_spec.set_position([.1, .1, .85, .6 ])
wv_spec = fig.axes[1]
wv_spec.set_position([.1, .8, .85, .15 ])
# freqs, power = signal.periodogram(data, fs=sample_rate, window='hanning', 
                                  # nfft=NFFT)
ax_spec.plot(freqs, power_db)
ax_spec.set_xscale('log')
ax_spec.set_xlabel('Frequency [Hz]')
ax_spec.set_ylabel('Power [dB ref. counts^2/Hz]')
ax_spec.set_ylim(-40, 60)
ax_spec.grid()

'''
# %%