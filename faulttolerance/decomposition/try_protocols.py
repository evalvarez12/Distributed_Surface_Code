import protocols
import numpy as np

epl = protocols.pair_EPL(ps=0,
                         pm=0,
                         pg=0,
                         eta=0,
                         a0=0,
                         a1=0,
                         theta=np.pi/4)
print(type(epl))
