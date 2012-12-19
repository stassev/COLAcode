CFLAGS= -std=c99 -O3 -Wall 
LIBS= -lm -lgsl -lgslcblas -lfftw3f 

COLAcode :  auxPM.o gadget_blob.o
		gcc ${CFLAGS} -o COLAcode main.c  auxPM.o gadget_blob.o ${LIBS}

auxPM.o :
		gcc ${CFLAGS} -c auxPM.c 

gadget_blob.o :
		gcc ${CFLAGS} -c gadget_blob.c 

clean:
	rm -f *.o COLAcode 
