/*******************************************************************************

 File:        os_al.c

 Description: This module provides an OS abstraction layer.
              It defines OS independant:
                   - Mutex (exclusive access to shared ressources)
				   - Semaphore (synchronization mechanism)
				   - Thread (creation, priority, etc...)
 
              This is the Windows NT implementation.
			  
 Author:      D. Calvet,        calvet@hep.saclay.cea.fr
              I. Mandjavidze, mandjavi@hep.saclay.cea.fr
			  J.M. Conte,     jmconte@hep.saclay.cea.fr

 History:
	April/98: function names changed

*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <semaphore.h>
#include <errno.h>
#include <limits.h>
#include <fcntl.h>
#include <signal.h>
#include <time.h> 
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/msg.h>
#include <sys/stat.h>

#include <string.h>

#include "Platform.h"
#include "Pipes.h"

/* Debug Section */
#define DBG_Pipe_DeleteInternal  
#define DBG_Pipe_CreateInternal  
#define DBG_Pipe_Create    
#define DBG_Pipe_Delete    
#define DBG_Pipe_Connect  
#define DBG_Pipe_Read  
#define DBG_Pipe_Write  

/*
 * Define Pipe structure
 */
typedef struct _Pipe
{
	int            msg_id;
} Pipe;
/* Types of nodes */
#define DEF_PIPE_MSG_MASTER 1
#define DEF_PIPE_MSG_SLAVE  0

/* Define max number of messages in queues and max size of messages */
#define DEF_PIPE_MAX_MSG_SIZE  4056


/******************************************************************************/
/* Pipe_Delete                                                                */
/*  0 on success                                                              */
/* -1 on failure                                                              */
/******************************************************************************/
int Pipe_DeleteInternal( void* *pi, int flag )
{
	Pipe *pipe;
	int status = 0;

	/* Chek for input parameters */
	if( !pi )
	{
		fprintf( stderr, "%s: wrong input parameter, pipe=0\n", __FUNCTION__ );
		return -1;
	}
	pipe = *pi;
	*pi = 0;
	if( !pipe )
	{
		fprintf( stderr, "%s: wrong input parameter, *pipe=0\n", __FUNCTION__ );
		return -1;
	}
	DBG_Pipe_DeleteInternal( "%s: pipe=0x%x\n", __FUNCTION__, pipe );

	/* Free message buffers */
	DBG_Pipe_DeleteInternal( "Pipe_DeleteInternal: flag=0x%x\n", flag );
	if( flag )
	{
		/* Free messages */
		if( pipe->msg_id != -1 )
		{
			if( msgctl( pipe->msg_id, IPC_RMID, (struct msqid_ds *)NULL ) == -1 )
			{
				fprintf( stderr, "%s: msgctl for msg_id failed\n", __FUNCTION__ );
				status = -1;
			}
			DBG_Pipe_DeleteInternal( "%s: del msg_id=%d\n", __FUNCTION__, pipe->msg_id );
		}
	}

	/* Free pipes structure */
	free( pipe );
	DBG_Pipe_DeleteInternal( "%s: pipe freed\n", __FUNCTION__ );

	return status;
}
int Pipe_Delete( void* *pi )
{
	/* Create pipes structure */
	if( Pipe_DeleteInternal( pi, 1 ) < 0 )
	{
		fprintf( stderr, "%s: Pipe_DeleteInternal failed\n", __FUNCTION__ );
		return -1;
	}
	return 0;
}
int Pipe_Close( void* *pi )
{
	/* Create pipes structure */
	if( Pipe_DeleteInternal( pi, 0 ) < 0 )
	{
		fprintf( stderr, "%s: Pipe_DeleteInternal failed\n", __FUNCTION__ );
		return -1;
	}
	return 0;
}

/******************************************************************************/
/* Pipe_Create                                                                */
/*  0 on success                                                              */
/* -1 on failure                                                              */
/******************************************************************************/
/*
 * Pipe_Create : Creates named pipe in server
 *  0 on success
 * -1 on failure
 */
int Pipe_CreateInternal( void* *pi, char *pipe_name, int flags )
{
	Pipe *pipe;
	key_t key;
	int msqflg; /* access permission */
	char filename[132];
	FILE *fptr;
	struct msqid_ds msg_ds;

	/* prepare output parameter for case of errors */
	*pi=0;

	/* Create pipes structure */
	if( ( pipe = (Pipe *)malloc( sizeof( Pipe ) ) ) == 0 )
	{
		fprintf( stderr, "%s: malloc failed for Pipe of %d bytes\n",
			__FUNCTION__, (int)sizeof( Pipe ) );
		return -1;
	}
	/* Clear the pipe structure */
	pipe->msg_id = -1;

	/* Create pipe structure */
	sprintf( filename, "/tmp/%s.msg", pipe_name );
	DBG_Pipe_CreateInternal( "%s: msg filename=%s\n",
		__FUNCTION__, filename );
	if( ( fptr = fopen( filename, "wr" ) ) == NULL )
	{
		fprintf( stderr, "%s: fopen failed for %s with %d %s\n",
			__FUNCTION__, filename, errno, strerror(errno) );
		Pipe_DeleteInternal( ( void* *)(&pipe), flags );
		return -1;
	}
	fclose( fptr );
	if( ( key = ftok( filename, 's' ) ) == -1 )
	{
		fprintf( stderr, "%s: ftok failed with %d %s\n", __FUNCTION__, errno, strerror(errno) );
		Pipe_DeleteInternal( ( void* *)(&pipe), flags );
		return -1;
	}
	DBG_Pipe_CreateInternal( "%s: key=%d\n", __FUNCTION__, key );
	msqflg = flags | 0666;
	if( (pipe->msg_id = msgget( key, msqflg ) ) == -1 )
	{
		fprintf( stderr, "%s: msqget for s2c_msg_id I failed with %d %d\n", __FUNCTION__, errno, strerror(errno) );
		msqflg =  0666;
		if( (pipe->msg_id = msgget( key, msqflg ) ) == -1 )
		{
			fprintf( stderr, "%s: msqget for msg_id II failed with %d %d\n", __FUNCTION__, errno, strerror(errno) );
			Pipe_DeleteInternal( ( void* *)(&pipe), flags );
			return( -1 );
		}
	}
	DBG_Pipe_CreateInternal( "%s: msg_id=%d\n", __FUNCTION__, pipe->msg_id );

	if( msgctl( pipe->msg_id, IPC_STAT, &msg_ds ) < 0 )
	{
		fprintf( stderr, "%s: msgctl IPC_STAT for msg_id failed with %d %s\n", __FUNCTION__, errno, strerror(errno) );
		Pipe_DeleteInternal( ( void* *)(&pipe), flags );
		return( -1 );
	}
	DBG_Pipe_CreateInternal( "%s: msg_qbytes=%d\n", __FUNCTION__, msg_ds.msg_qbytes );

	/* set output parameter */
	*pi = pipe;

	return ((Pipe *)*pi)->msg_id;
}
int Pipe_Create( void* *pi, char *pipe_name, int direction )
{
	int ret;
	/* Create pipes structure */
	if( (ret = Pipe_CreateInternal( pi, pipe_name, IPC_CREAT )) < 0 )
	{
		fprintf( stderr, "%s: Pipe_CreateInternal failed\n", __FUNCTION__ );
		return -1;
	}
	return ret;
}

int Pipe_Open( void* *pi, char *pipe_name, int direction )
{
	int ret;
	/* Create pipes structure */
	if( (ret = Pipe_CreateInternal( pi, pipe_name, 0 )) < 0 )
	{
		fprintf( stderr, "%s: Pipe_CreateInternal failed\n", __FUNCTION__ );
		return -1;
	}
	return ret;
}


/******************************************************************************/
/* Pipe_Connect                                                               */
/*  0 on success                                                              */
/* -1 on failure                                                              */
/******************************************************************************/
int Pipe_Connect( void *pi )
{
	return 0;
}


/******************************************************************************/
/* Pipe_Read                                                                  */
/*  0 on success                                                              */
/* -1 on failure                                                              */
/******************************************************************************/
int Pipe_Read( void *pi, char *buf, int *sz)
{
	Pipe *pipe;
	struct msgbuf *msg_buf;
	int rcv;

	/* Chek for input parameters */
	pipe = (Pipe *)pi;
	if( !pipe )
	{
		fprintf( stderr, "%s: wrong input parameter, *pipe=0\n", __FUNCTION__ );
		return -1;
	}
	DBG_Pipe_Read( "%s: pipe=0x%x buf=0x%x sz=0x%x\n", __FUNCTION__, pipe, buf, sz );

	/* Receive data */
	msg_buf = (struct msgbuf *)buf;
	if( (rcv = msgrcv( pipe->msg_id, msg_buf, DEF_PIPE_MAX_MSG_SIZE, 0, MSG_NOERROR ) )
	< 0 )
	{
		fprintf( stderr, "%s: msgrcv faild for msg_id %d failed with %d %s\n", __FUNCTION__, pipe->msg_id, errno, strerror(errno) );
		return( -1 );
	}

	/* Copy received data to output */
	*sz = rcv;
	DBG_Pipe_Read( "%s: buffer of %d bytes received\n", __FUNCTION__, *sz );
	return( 0 );
}

/******************************************************************************/
/* Pipe_Write                                                                 */
/*  0 on success                                                              */
/* -1 on failure                                                              */
/******************************************************************************/
int Pipe_Write( void *pi, char *buf, int sz)
{
	Pipe *pipe;
	struct msgbuf *msg_buf;
	int snd;
	/* Chek for input parameters */
	pipe = (Pipe *)pi;
	if( !pipe )
	{
		fprintf( stderr, "%s: wrong input parameter, *pipe=0\n", __FUNCTION__ );
		return -1;
	}
	DBG_Pipe_Write( "%s: pipe=0x%x buf=0x%x sz=0x%x\n", __FUNCTION__, pipe, buf, sz );
	if( sz > DEF_PIPE_MAX_MSG_SIZE )
	{
		fprintf( stderr, "%s: message sz=%d > MAX_MSG_SIZE=%d\n",
			__FUNCTION__, sz, DEF_PIPE_MAX_MSG_SIZE );
		return -1;
	}

	/* Send data */
	msg_buf = (struct msgbuf *)buf;
	if( msgsnd( pipe->msg_id, msg_buf, sz, 0 ) < 0 )
	{
		fprintf( stderr, "%s: msgsnd faild for msg_id %d failed with %d %s\n", __FUNCTION__, pipe->msg_id, errno, strerror(errno) );
		return( -1 );
	}
	DBG_Pipe_Write( "%d: Done\n", __FUNCTION__ );
	return 0;
}
