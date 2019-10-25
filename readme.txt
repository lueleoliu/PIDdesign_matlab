Station.txt格式:
is_pitch	Station	is_first	first_td

config.txt格式：
Generations	PoulationSize	Fit_Coefficient

Result.txt格式：
for Pitch：
Kp	Ti	Kd	Td（Same as Pitch1）	SettlingTime	RiseTime		Overshoot	Oscnum
for Torque：
Kp	Ti	SettlingTime	RiseTime

baseline.txt格式：
Pm_base		Wpm_base	Gm_base

range.txt格式：
Pitch_Kp_low_band		Pitch_Ti_low_band		Pitch_Kd_low_band		Pitch_Ti_low_band
Pitch_Kp_high_band		Pitch_Ti_high_band		Pitch_Kd_high_band		Pitch_Ti_high_band
Torque_Kp_low_band	Torque_Ti_low_band	0			0
Torque_Kp_high_band	Torque_Ti_high_band	0			0
