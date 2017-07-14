import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy.signal import butter, lfilter, freqz

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float) ** 2
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def butter_highpass(cutoff, fs, order=3):
    nyq = 0.5 * fs
    b, a = butter(order, cutoff / nyq, btype='high')
    return b, a

def highpass(data, cutoff, fs, order=3):
    b, a = butter_highpass(cutoff, fs, order)
    return lfilter(b, a, data)

def butter_bandpass(lowcut, highcut, fs, order=3):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def bandpass(data, lowcut, highcut, fs, order=3):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    return lfilter(b, a, data)

def band_powers(data, fs, freqs, pulse_duration,
        sample_duration, sweep_duration, find_bandwidth=400, total_bandwidth=5000, debug=False):
    channels = xrange(data.shape[1])
    pulse_length = int(fs * pulse_duration)
    sample_length = int(fs * sample_duration)
    sweep_length = int(fs * sweep_duration)
    tones = []
    powers = {}
    for channel in channels:
        data[:,channel] = highpass(data[:,channel], 20e3, fs, 5)
    # Find and remove the frequency sweep
    sweep_index = len(data)
    for channel in channels:
        sweep_index = min(
            sweep_index,
            np.argmax(moving_average(data[:,channel], sweep_length))
        )
    data = data[0:sweep_index]
    for freq in freqs:
        avg_powers = []
        for channel in channels:
            band = bandpass(data[:,channel], freq - find_bandwidth / 2, freq + find_bandwidth / 2, fs)
            tone_start = np.argmax(moving_average(band, pulse_length))
            if channel == 0:
                tones.append(tone_start)
            band = bandpass(data[:,channel], freq - total_bandwidth / 2, freq + total_bandwidth / 2, fs)
            power = np.sum(band[tone_start:tone_start+sample_length] ** 2)
            print "POWER", power
            avg_powers.append(power)
        powers[freq] = np.mean(avg_powers)
    if debug:
        plt.plot(np.array(tones, dtype=np.float) / fs, freqs, 'bx')
        plt.specgram(data[:,0], Fs=fs, NFFT=1024)
        plt.show()
    print "POWERS", powers
    return powers