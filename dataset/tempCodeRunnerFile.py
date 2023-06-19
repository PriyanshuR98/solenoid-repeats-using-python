import numpy as np
import matplotlib.pyplot as plt


dicfact1={}
dicfact1={'A':-0.591 , 'C':-1.343,'D':1.050 , 'E':1.375 , 'F':-1.006 , 'G':-0.384 , 'H':0.036 ,'I':-1.29,'K':1.831,'L':-1.019,'M':-0.663,'N':0.945,'P':0.189,'Q':0.931,'R':1.538,'S':-0.228,'T':-0.032,'V':-1.337,'W':-0.595,'Y':0.260}

dicfact2={}
dicfact2={'A':-1.302,'C':0.465,'D':0.302 , 'E':-1.453 , 'F':-0.590 , 'G':1.652 , 'H':-0.417 ,'I':-0.547,'K':-0.561,'L':-0.987,'M':-1.524,'N':0.828,'P':2.081,'Q':-0.179,'R':-0.055,'S':1.399,'T':0.326,'V':-0.279,'W':0.009,'Y':0.830}

dicfact3={}
dicfact3={'A':-0.733,'C':-0.862,'D':-3.656 , 'E':1.477 , 'F':1.891 , 'G':1.330 , 'H':-1.673 ,'I':2.131,'K':0.533,'L':-1.505,'M':2.219,'N':1.299,'P':-1.628,'Q':-3.005,'R':1.502,'S':-4.760,'T':2.213,'V':-0.544,'W':0.672,'Y':3.097}

dicfact4={}
dicfact4={'A':1.570,'C':-1.020,'D':-0.259 , 'E':0.113 , 'F':-0.397 , 'G':1.045 , 'H':-1.474 ,'I':0.393,'K':-0.277,'L':1.266,'M':-1.005,'N':-0.169,'P':0.421,'Q':-0.503,'R':0.440,'S':0.670,'T':0.908,'V':1.242,'W':-2.128,'Y':-0.838}

dicfact5={}
dicfact5={'A':-0.146,'C':-0.255,'D':-3.242 , 'E':-0.837 , 'F':0.412 , 'G':2.064 , 'H':-0.078 ,'I':0.816,'K':1.648,'L':-0.912,'M':1.212,'N':0.933,'P':-1.392,'Q':-1.853,'R':2.897,'S':-2.647,'T':1.313,'V':-1.262,'W':-0.184,'Y':1.512}

s="GAMGISLGNSEADRQLLEAAKAGDVETVKKLCTVQSVNCRDIEGRQSTPLHFAAGYNRVSVVEYLLQHGADVHAKDKGGLVPLHNACSYGHYEVAELLVKHGAVVNVADLWKFTPLHEAAAKGKYEICKLLLQHGADPTKKNRDGNTPLDLVKDGDTDIQDLLRGDAALLDAAKK"
# print(len(s))
ans1=""
vec1=list();
for i in range(len(s)):
    ans1+=str(dicfact1[s[i]])+" "
    vec1.append((dicfact1[s[i]]))
    
ans2=""
vec2=list();
for i in range(len(s)):
    ans2+=str(dicfact2[s[i]])+" "
    vec2.append((dicfact2[s[i]]))

ans3=""
vec3=list();
for i in range(len(s)):
    ans3+=str(dicfact3[s[i]])+" "
    vec3.append((dicfact3[s[i]]))
    
ans4=""
vec4=list();
for i in range(len(s)):
    ans4+=str(dicfact4[s[i]])+" "
    vec4.append((dicfact4[s[i]]))
    
ans5=""
vec5=list();
for i in range(len(s)):
    ans5+=str(dicfact5[s[i]])+" "
    vec5.append((dicfact5[s[i]]))

x1=np.array(vec1)
op1=np.fft.fft(x1)
N1=len(op1)
freq1=np.arange(N1)
# plt.plot(freq, np.abs(op))
# plt.show()
x2=np.array(vec2)
op2=np.fft.fft(x2)
N2=len(op2)
freq2=np.arange(N2)
#plt.plot(freq, np.abs(op))
#plt.show()
x3=np.array(vec3)
op3=np.fft.fft(x3)
N3=len(op3)
freq3=np.arange(N3)
#plt.plot(freq, np.abs(op))
#plt.show()
x4=np.array(vec4)
op4=np.fft.fft(x4)
N4=len(op4)
freq4=np.arange(N4)
#plt.plot(freq, np.abs(op))
#plt.show()
x5=np.array(vec5)
op5=np.fft.fft(x5)
N5=len(op5)
freq5=np.arange(N5)
#plt.plot(freq, np.abs(op))
#plt.show()

fig,axis_plot=plt.subplots(nrows=5, ncols=0, figsize=(10, 12))
axis_plot[0][0].plot(freq1, np.abs(op1))
axis_plot[0][0].set_title('Factor 1')
axis_plot[1][0].plot(freq2, np.abs(op2))
axis_plot[1][0].set_title('Factor 2')
axis_plot[2][0].plot(freq3, np.abs(op3))
axis_plot[2][0].set_title('Factor 3')
axis_plot[3][0].plot(freq4, np.abs(op4))
axis_plot[3][0].set_title('Factor 4')
axis_plot[4][0].plot(freq5, np.abs(op5))
axis_plot[4][0].set_title('Factor 5')


fig.suptitle('Protein Sequence DFT')


plt.show()