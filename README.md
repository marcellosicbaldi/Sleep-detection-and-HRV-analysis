# Sleep-detection-and-HRV-analysis
Implementation of a sleep detection algorithm using accelerometer data, followed by an evaluation of HRV features during sleep and wake using PPG data. Both PPG and accelerometer data were recorded in 10 healthy subjects using the Empatica E4 wrist-wearable device.

- [project report](https://github.com/marcellosicbaldi/Sleep-detection-and-HRV-analysis/blob/main/project_report_Marcello_Sicbaldi.pdf) 

## Sleep detection

Periods of sleep and wake are detected using an algorithm based on thresholding the Z-angle computed from accelerometer data. Here is an example for one subject:
![step8](https://user-images.githubusercontent.com/55695116/174432828-2203c5ac-8ed1-42b7-bb7a-834f69acb3a9.jpg)

## HRV comparison between sleep and wake
 
The analysis is based on inter-beat intervals computed from the PPG signal. The HRV features I decided to investigate are the Heart Rate (HR), standard deviation of normal beats (SDNN), root mean square of successive differences (RMSSD), low and high frequency power and their ratio. 

Here are some of the results:

![pretty_plot](https://user-images.githubusercontent.com/55695116/174432973-3c0fab41-cfce-49f7-80dc-5f5711597374.png)
![RMSSD](https://user-images.githubusercontent.com/55695116/174432976-ccd61d39-19f5-4d7f-9572-332ede347a64.png)
![SDNN](https://user-images.githubusercontent.com/55695116/174432978-dd8e74f0-edc7-46ce-8049-31a3659adbaf.png)


