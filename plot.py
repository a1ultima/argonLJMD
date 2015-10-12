#!/usr/bin/env python

import matplotlib.pyplot as plt

fi = open("./data/temperature_diffusion-coefficient.dat","r")

dc_list = []
dc2_list = []
dc3_list = []
dc4_list = []
dt_c_list = []

dT_list = []
dt_list = []


while True:

	line = fi.readline()

	if line == "":
		break

	if "#" in line:
		continue

	line_split = line.split("\t")

	dt = float(line_split[0])

	dc = float(line_split[2])
	dc_list.append(dc)
	
	dT = float(line_split[1])
	dt_list.append(dt)
	dT_list.append(dT)


# red dashes, blue squares and green triangles
plt.plot(dt_list, dc_list)

fig = plt.gcf()

plt.ylabel('Diffusion Coefficient')
plt.xlabel('Time')
plt.title('Diffusion Coefficient of Argon particles')
plt.savefig('./plots/diffusion_coefficient.png')

plt.plot(dt_list[1:], dT_list[1:], 'g')

fig = plt.gcf()

plt.ylabel('Temperature')
plt.xlabel('Time')
plt.title('Temperature fluctuations of Argon particles')
plt.savefig('./plots/temperature.png')

fi.close()

