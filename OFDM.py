
import commpy 
from commpy.modulation import PSKModem, QAMModem
import numpy as np

def modulate(BasebandSignal,type):
    if(type=='4PSK'):
        modulator=PSKModem(4)
    elif (type=='8PSK'):
        modulator=PSKModem(8)
    elif (type=='16PSK'):
        modulator=PSKModem(16)
    elif (type=='4QAM'):
        modulator=QAMModem(4)
    elif (type=='8QAM'):
        modulator=QAMModem(8)
    else:
        modulator=QAMModem(16)
    return modulator,modulator.modulate(BasebandSignal)


def sample(Signal,k,carr_num,sign):
    x_real=0
    x_img=0
    for i in range(carr_num):
        x_real +=Signal[i].real*np.cos(2*np.pi*i*k/carr_num)-Signal[i].imag*np.sin(2*sign*np.pi*i*k/carr_num)
        x_img +=Signal[i].real*np.sin(2*sign*np.pi*i*k/carr_num)+Signal[i].imag*np.cos(2*np.pi*i*k/carr_num)
    return (complex(x_real,x_img))/np.sqrt(carr_num)


def OFDM_encode(Signal,H,Symbol_rate,Bandwidth):
    # Assume the cyclic prefix is equal to maximum delay spread
    max_delay=len(H)-1
    #Calculate the number of sub-carrier
    Carr_num=np.ceil(max_delay/(Bandwidth/Symbol_rate-1))
    Carr_num=int(Carr_num)
    # Do sampling in every window
    block_index=0
    block_num=int(len(Signal)/Carr_num)
    encoded_sig=np.zeros((block_num,Carr_num+len(H)-1),dtype=complex)
    trans_sig=np.zeros((block_num,Carr_num),dtype=complex)
    # Do OFDM encoder
    while block_index<(block_num):
        # Current block and next block
        s_cu=Signal[block_index*Carr_num:(block_index+1)*Carr_num]
        s_next=Signal[(block_index+1)*Carr_num:(block_index+2)*Carr_num]
        for t in range(Carr_num):
            # Do sampling
            encoded_sig[block_index][t+len(H)-1]=sample(Signal=s_cu,k=t,carr_num=Carr_num,sign=1)
            # Overlap
            # Discard surplus data
            if (t>Carr_num-len(H) and len(s_next)==Carr_num):
                encoded_sig[block_index+1][t-Carr_num+len(H)-1]=sample(Signal=s_next,k=(t-Carr_num),carr_num=Carr_num,sign=1)
        block_index +=1
    # encoded signal cross a multipath channel
    for i in range(block_num):
        # ignore cyclic prefix
        for j in range(Carr_num):
            for k in range(len(H)):
                trans_sig[i][j] +=H[k]*encoded_sig[i][j-k+len(H)-1]
            j=j+1
    
    return Carr_num,trans_sig

def OFDM_Decode(trans_sig,H,carr_num):
    channel_index=[]
    for i in range(carr_num):
        if(i<len(H)):
            channel_index.append(H[i])
        else:
            channel_index.append(0)
    channel_index=np.fft.fft(channel_index)
    decode_sig=[]
    for t in range(np.shape(trans_sig)[0]):
        for k in range(carr_num):
            temp=sample(Signal=trans_sig[t],k=k,carr_num=carr_num,sign=-1)
            temp=temp/channel_index[k]
            decode_sig.append(temp)
    print(len(decode_sig))
    return decode_sig


    

if __name__=='__main__':
    BasebandSignal=[]
    Seq_len=76
    for i in range(Seq_len):
        BasebandSignal.append(np.random.randint(0,2))
    modulator,BandpassSignal=modulate(BasebandSignal,'8PSK')
    # Parameter of channel
    H=[0.5,0.6]
    # Desired Symbol rate 
    Symbol_rate=5e+05
    #available bandwith
    Bandwidth=6e+05
    carrnum,tran_sig=OFDM_encode(BandpassSignal,H,Symbol_rate,Bandwidth)
    decode_sig=OFDM_Decode(tran_sig,H,carrnum)
    demodulate_sig=modulator.demodulate(decode_sig,'hard')