/*******************************************************************************

 File:        synchro.h

 Description: This module provides an OS abstraction layer.
              It declares OS independant:
		- Mutex (exclusive access to shared ressources)
		- Semaphore (synchronization mechanism)
		- Thread (creation, priority, etc...)
              See synchro.c in the OS specific directory for implementation.

 Author:      D. Calvet,        calvet@hep.saclay.cea.fr
              I. Mandjavidze, mandjavi@hep.saclay.cea.fr

 History:
 	Apr/98: function names changed
 	   
*******************************************************************************/

#ifndef PIPE_H
#define PIPE_H

/* Operating system independent bi-directionnal named pipes for messages */
#define DEF_Pipe_Type_Wr 0
#define DEF_Pipe_Type_Rd 1
int Pipe_Create(  void* *pi, char *pipe_name, int direction );
int Pipe_Delete(  void* *pi );
int Pipe_Connect( void  *pi );
int Pipe_Open(    void* *pi, char *pipe_name, int direction );
int Pipe_Close(   void* *pi );
int Pipe_Read(    void  *pi, char *buf, int *sz);
int Pipe_Write(   void  *pi, char *buf, int  sz);

#endif /* #ifndef PIPE_H */
