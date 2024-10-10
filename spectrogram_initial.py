
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

data_path = '/data/stor/basic_data/seismic_data/day_vols/WOLVERINE/WOLC/'
file_prefix = 'WOLC.XX.00.HHZ.2022.' 
day = 244 #122 152 182 213 244 
st = obspy.read(data_path + file_prefix + str(day), format='MSEED')



#%%
pre_filt1 = [0.05, 0.1, 80, 100] #standard bandpass parameters to remove low and high freq. signal 
instr_resp = read_inventory("/data/stor/basic_data/seismic_data/day_vols/WOLVERINE/resp/Wolverine_station_20240824.xml")
#%%
start = st[0].stats.starttime
end = st[0].stats.endtime
station_cut = st.copy()
station_cut.trim(starttime=start+(20*60*60), endtime = start+(20*60*60)+60*60)
station_rem_spec = station_cut.copy()
station_rem_spec.remove_response(inventory=instr_resp, water_level= None,pre_filt=pre_filt1)

#%% Spectrogram with waveform too.
fig, ax = plt.subplots(num=3, clear=True, figsize=(16,10))
station_rem_spec.plot(fig=fig) # Plot waveform into the figure fig

# Define and reshape axes
ax_spec = fig.axes[0] 
ax_spec.set_position([.1, .1, .85, .6 ])

wv_spec = fig.axes[1]
wv_spec.set_position([.1, .8, .85, .15 ])
#wv_spec.set_ylim(-15000, 15000)
# wv_spec.set_ylim(-4000, 4000) # change these waveform axes limits to suit!

# freqs, power = signal.periodogram(data, fs=sample_rate, window='hanning', 
                                  # nfft=NFFT)
# NFFT is the number of data points to use in calculating periodograms.
#   It's usually a power of 2.  Smaller numbers provide better temporal resolution
#   and larger numbers provide better (low frequency) spectral resolution.
NFFT = 2**2
NFFT = 2**14
# Calculate the spectrogram.  vmin and vmax define the power bounds of the colormap
ax_spec.specgram(station_rem_spec[0].data, NFFT=NFFT, Fs=int(station_rem_spec[0].stats.sampling_rate), Fc=None, 
                 detrend='linear', window=None, noverlap=NFFT*.5, cmap=None, 
                 xextent=None, pad_to=None, sides=None, scale_by_freq=True, 
                 mode='psd', scale=None)#, vmin=None, vmax=-40)  # *** how do i define v min v max 
                #  mode='psd', scale=None, vmin=-15, vmax=45)
                 #mode='psd', scale=None, vmin=-20, vmax=30)
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
# %%
