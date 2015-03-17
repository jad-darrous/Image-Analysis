#include <stdlib.h>
#include <stdio.h>
#include "utils.h"
#include "pixmapIO.h"

int * readPixmap(char* filename,
			   int * type,
			   int * nbColumns,
			   int * nbRows)
{
  FILE * filePtr;
  int* pixmap;
  int ich1, ich2;
  int rows, cols ;
  int maxVal;
  int i,j,k;

  /* Additional varaibles */
  int read8Bit;
  int singleBit;
  int mask;

  
  /* File opening */
  filePtr = fopen(filename,"rb");
  
  /* Magic number reading */
  ich1 = getc( filePtr );
  if ( ich1 == EOF )
    pm_erreur( "EOF / read error reading magic number" );
  ich2 = getc( filePtr );
  if ( ich2 == EOF )
    pm_erreur( "EOF / read error reading magic number" );

  switch(ich2) 
    {
    default:
      fclose(filePtr);
      printf("File format not supported\n");
      exit(1);
      break;
      
    case '1':
      *type = 1;

      /* Reading dimensions */
      cols = pm_getint( filePtr );
      rows = pm_getint( filePtr );
      /* Passing number of columns and number of rows as arguments */
      *nbColumns = cols;
      *nbRows = rows;
      
      /* Memory allocation */
      pixmap = (int *) malloc(cols * rows * sizeof(int));
      
      /* Reading data */
      for(i=0; i < rows; i++)
	for(j=0; j < cols ; j++)
	  pixmap[i * cols + j] = pm_getbit(filePtr);
      break;

    case '2':
      *type = 2;
      /* Reading dimensions */
      cols = pm_getint( filePtr );
      rows = pm_getint( filePtr );
      maxVal = pm_getint( filePtr );
      /* Passing number of columns and number of rows as arguments */
      *nbColumns = cols;
      *nbRows = rows;
      
      /* Memory allocation */
      pixmap = (int *) malloc(cols * rows * sizeof(int));
      
      /* Reading data */
      for(i=0; i < rows; i++)
	for(j=0; j < cols ; j++)
	  pixmap[i * cols + j] = (int) pm_getint(filePtr);
      break;


    case '3':
      *type = 3;
     
      /* Reading dimensions */
      cols = pm_getint( filePtr );
      rows = pm_getint( filePtr );
      maxVal = pm_getint( filePtr );
      /* Passing number of columns and number of rows as arguments */
      *nbColumns = cols;
      *nbRows = rows;
      printf("%d %d \n", cols, rows);
      printf("%d \n", maxVal);
      /* Memory allocation */
      pixmap = (int *) malloc(cols * rows * 3 * sizeof(int));
      
      /* Reading data */
  
      for(i=0; i < rows; i++)
	for(j=0; j < cols ; j++)
	  for(k=0 ; k < 3 ; k++)  
	    pixmap[i * cols * 3 + j * 3 + k] = (int) pm_getint(filePtr);

      break;




    case '4':
      *type = 4;
      printf("Type: 4\n");
      /* Reading dimensions */ 
      cols = pm_getint( filePtr );
      rows = pm_getint( filePtr );

      /* Passing number of columns and number of rows as arguments */
      *nbColumns = cols;
      *nbRows = rows;
      printf("cols: %d, rows: %d\n", cols, rows);

      /* Memory allocation */
      pixmap = (int *) malloc(rows * cols * sizeof(int));
      
      /* Reading data  bit by bit... tough !*/
      /* http://orion.math.iastate.edu/burkardt/g_src/pbmpak/pbmpak.c */

      for(j=0; j < rows; j++) {
	for(i=0; i < cols ; i++) {
	  if (i%8 == 0) {
	    read8Bit = fgetc(filePtr);
	    
	    if (read8Bit == EOF) 
	      pm_erreur("EOF / read error, file dimension is not a power of 2");

	  }

	  mask = 7 - i%8;
	  singleBit = ( read8Bit >> mask )%2;

	  pixmap[j * cols + i] = (int) singleBit;
	}
      }
      
      break;

    case '5':
      *type = 5;
      /* Reading dimensions */
      cols = pm_getint( filePtr );
      rows = pm_getint( filePtr );
      maxVal = pm_getint( filePtr );
      /* Passing number of columns and number of rows as arguments */
      *nbColumns = cols;
      *nbRows = rows;
      
      /* Memory allocation */
      pixmap = (int *) malloc(cols * rows * sizeof(int));
      
      /* Reading data */
      for(i=0; i < rows; i++)
	for(j=0; j < cols ; j++)
	  pixmap[i * cols + j] = pm_getrawbyte(filePtr);

      break;

    case '6':
      *type = 6;
      /* Reading dimensions */
      cols = pm_getint( filePtr );
      rows = pm_getint( filePtr );
      maxVal = pm_getint( filePtr );
      /* Passing number of columns and number of rows as arguments */
      *nbColumns = cols;
      *nbRows = rows;
      
      /* Memory allocation */
      pixmap = (int *) malloc(cols * rows * 3 * sizeof(int));
      
      /* Reading data */
      for(i=0; i < rows; i++)
	for(j=0; j < cols ; j++)
	  for(k=0 ; k < 3 ; k++)
	    pixmap[i * cols * 3 + j * 3 + k] =  pm_getrawbyte(filePtr);
      break;
      
    }
  
  fclose(filePtr);
  return pixmap;
}




int writePixmap(int * pixmap,
		 int cols, int rows,
		 int type,
		 char * filename)
{
  int i, j, k;
  FILE * filePtr;

  /* Additional variables for binary file writing*/
  int bit;
  unsigned char c;
  int blackOrWhite;
  int mask;

  /* Provides the right extension */
  char * extension =  (char *) malloc(5*sizeof(char));
  if ( (type == 1) || (type == 4))
    extension = "pbm";
  else if ( (type == 2) || (type == 5) )
    extension = "pgm";
  else if ( (type == 3) || (type == 6) )
    extension = "ppm";
  else
    extension = 0;
  pm_setExtension(filename, extension);


  filePtr = fopen(filename, "wb");

  switch (type)
    {

    case 1:
      fprintf(filePtr, "P1\n");
      fprintf(filePtr, "%d %d \n", cols, rows);
      
      for (j = 0; j < rows; j++) {
	for (i = 0; i < cols; i++) {
	  if (pixmap[j*cols + i] > 0)
	    fprintf(filePtr, "1\n");
	  else
	    fprintf(filePtr, "0\n");
	}
      }  
      fprintf(filePtr, "\n");

      break;

    case 2:
      fprintf(filePtr, "P2\n");
      fprintf(filePtr, "%d %d \n", cols, rows);
      fprintf(filePtr, "%u\n", (unsigned char) pm_getMax(pixmap, cols * rows));
      
      for (j = 0; j < rows; j++)
	for (i = 0; i < cols; i++)
	  fprintf(filePtr, "%u\n", (unsigned char) pixmap[j*cols + i]);
      fprintf(filePtr, "\n");

      break;


      case 3:

	fprintf(filePtr, "P3\n");
	fprintf(filePtr, "%d %d \n", cols, rows);
	
	fprintf(filePtr, "%u\n", (unsigned char) pm_getMax(pixmap, cols * rows * 3));
	for (j = 0; j < rows; j++)
	  for (i = 0; i < cols; i++)
	    for (k = 0 ; k != 3 ; k++)
	      fprintf(filePtr, "%u\n", (unsigned char) pixmap[j*cols * 3+ i*3 + k]);
	fprintf(filePtr, "\n");
	
      break;
      


    case 4:
      /* Tough ! see http://orion.math.iastate.edu/burkardt/g_src/pbmpak/pbmpak.c */
      fprintf(filePtr, "P4\n");
      fprintf(filePtr, "%d %d \n", cols, rows);

      c = 0;
      for (j = 0; j < rows; j++) {
	for (i = 0; i < cols; i++) {
	  if (pixmap[j*cols + i] > 0)
	    blackOrWhite = 1;
	  else 
	    blackOrWhite = 0;

	  mask = 7 - i%8;
	  
	  bit = blackOrWhite % 2;
	  c = c | (bit << mask);
	  
	  if ( (i+1)%8 == 0 || i == (cols - 1) ) {
	    fputc(c, filePtr);
	    c = 0;
	  }
	  
	}
      }
      fprintf(filePtr, "\n");
      
      break;
      
    case 5:
    printf("-----\n");
      fprintf(filePtr, "P5\n");
      fprintf(filePtr, "%d %d \n", cols, rows);
      fprintf(filePtr, "%u\n", (unsigned char) pm_getMax(pixmap, cols * rows));
      
      for (j = 0; j < rows; j++)
	for (i = 0; i < cols; i++)
	  fputc((unsigned char) pixmap[j*cols + i], filePtr);
      fprintf(filePtr, "\n");

      break;

      
    case 6:
      fprintf(filePtr, "P6\n");
      fprintf(filePtr, "%d %d \n", cols, rows);
      fprintf(filePtr, "%u\n", (unsigned char) pm_getMax(pixmap, cols * rows *3));
      
      for (j = 0; j < rows; j++)
	for (i = 0; i < cols; i++)
	  for(k=0; k !=3 ; k++)
	    fputc((unsigned char) pixmap[j*cols*3 + i*3 +k], filePtr);
      fprintf(filePtr, "\n");

      break;

    default:
      fclose(filePtr);
      pm_erreur("File type not supported yet...\n");
      break;
    }

  fclose(filePtr);
  return 0;
}
