"""Generate a literal straight-line bsc-ops-v1 instruction list.

The Python generator/interpreter is not part of the RAM program. Each emitted
instruction is exactly one listed RAM primitive (except free register moves).
All schedule loops are expanded before the received word is examined.
"""
from construction import T48,T32,LUT48,LUT32,END

# A register operand is a one-tuple. Integer operands are literal constants.
def program():
    p=[]
    def emit(op,dst,*args):p.append((op,dst,args))
    def reg(i):return (i,)
    for block in range(42):
        n,k,start,dest,table,lut = ((48,24,48*block,24*block,T48,LUT48)
            if block<41 else (32,16,1968,984,T32,LUT32))
        word,off=divmod(start,64)
        emit('load',0,word)
        if off:emit('shr',0,reg(0),off)
        if off+n>64:
            emit('load',1,word+1);emit('shl',1,reg(1),64-off)
            emit('or_',0,reg(0),reg(1))
        emit('and_',0,reg(0),(1<<n)-1)
        emit('move',2,0)
        for c in range(n//8):
            emit('shr',1,reg(0),8*c);emit('and_',1,reg(1),255)
            emit('add',1,lut+c*256,reg(1));emit('load',1,reg(1))
            emit('xor',2,reg(2),reg(1))
        emit('add',2,table,reg(2));emit('load',2,reg(2))
        emit('xor',0,reg(0),reg(2));emit('and_',0,reg(0),(1<<k)-1)
        ow,oo=divmod(dest,64)
        emit('load',1,32+ow);emit('shl',2,reg(0),oo)
        emit('or_',1,reg(1),reg(2));emit('store',None,32+ow,reg(1))
        if oo+k>64:
            emit('load',1,33+ow);emit('shr',2,reg(0),64-oo)
            emit('or_',1,reg(1),reg(2));emit('store',None,33+ow,reg(1))
    return p

PROGRAM=program()
W_ACTUAL=sum(op!='move' for op,_,_ in PROGRAM)

def decode(machine):
    regs=[0,0,0]
    for op,dst,operands in PROGRAM:
        args=[regs[a[0]] if isinstance(a,tuple) else a for a in operands]
        if op=='move':regs[dst]=args[0]
        elif op=='store':machine.store(*args)
        else:regs[dst]=getattr(machine,op)(*args)

class TableMemory(dict):
    """Compact backing of fully preloaded fixed RAM words, not a new RAM op.

    Input/output words use the dictionary. All fixed table values are already
    computed in arrays before the RAM decoder starts. Machine.load is unchanged.
    """
    def __init__(self,initial,c48,c32,l48,l32):
        super().__init__(initial)
        self.c48,self.c32,self.l48,self.l32=c48,c32,l48,l32
    def get(self,address,default=0):
        if T48<=address<T32:return self.c48[address-T48]
        if T32<=address<LUT48:return self.c32[address-T32]
        if LUT48<=address<LUT32:return self.l48[address-LUT48]
        if LUT32<=address<END:return self.l32[address-LUT32]
        return super().get(address,default)
