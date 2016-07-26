IP=192.168.1.32

sshpass -p 'raspberry' scp -r pi@$IP:heat/sensors/rpi/temps .
