#include<stdio.h>
#include<stdlib.h>
#include<string.h>



int numberOfCells;
#define CELLS 1000000
int cells[CELLS];
int cTypes[CELLS];
int cStatus[CELLS];
void
read ()
{
	numberOfCells = 11659;
	int i, j;
	for (i = 0; i < numberOfCells; i++)
	  {
		  for (j = 0; j < 3; j++)
			  scanf ("%d ", &cells[i * 3 + j]);
		  scanf ("%d ", &cStatus[i]);

	  }
	for (i = 0; i < numberOfCells; i++)
		cTypes[i] = 1;
}


void
swap_Nbyte (char *data, int n, int m)
{
	int i, j;
	char old_data[16];

	for (j = 0; j < n; j++)
	  {
		  memcpy (&old_data[0], &data[j * m], m);
		  for (i = 0; i < m; i++)
			  data[j * m + i] = old_data[m - i - 1];
	  }
}

void
writeVTK ()
{
	FILE *fhandle;
	char header[256];
	int i;

	fhandle = fopen ("result.vtk", "wb");

	sprintf (header,
		 "# vtk DataFile Version 2.0\nStemscan output\nBINARY\nDATASET UNSTRUCTURED_GRID\n");
	fwrite (header, sizeof (char), strlen (header), fhandle);
	memset (header, 0, 256);
	sprintf (header, "\nPOINTS %d int\n", numberOfCells);
	fwrite (header, sizeof (char), strlen (header), fhandle);
	swap_Nbyte ((char *) cells, numberOfCells * 3, sizeof (int));
	fwrite (cells, sizeof (int), numberOfCells * 3, fhandle);
	memset (header, 0, 256);
	sprintf (header, "\nCELL_TYPES %d\n", numberOfCells);
	fwrite (header, sizeof (char), strlen (header), fhandle);
	swap_Nbyte ((char *) cTypes, numberOfCells, sizeof (int));
	fwrite (cTypes, sizeof (int), numberOfCells, fhandle);
	memset (header, 0, 256);
	sprintf (header, "\nPOINT_DATA %d", numberOfCells);
	fwrite (header, sizeof (char), strlen (header), fhandle);
	memset (header, 0, 256);
	sprintf (header,
		 "\nSCALARS status int 1\nLOOKUP_TABLE default\n");
	fwrite (header, sizeof (char), strlen (header), fhandle);
	swap_Nbyte ((char *) cStatus, numberOfCells, sizeof (int));
	fwrite (cStatus, sizeof (int), numberOfCells, fhandle);
	memset (header, 0, 256);
	//sprintf(header,"\nSCALARS size int 1\nLOOKUP_TABLE default\n");
	//fwrite(header,sizeof(char),strlen(header),fhandle);

	/*swap_Nbyte((char*)cSize,numberOfCells,sizeof(int));
	   fwrite(cSize,sizeof(int),numberOfCells,fhandle); */

	fclose (fhandle);
}

int
main ()
{
	read ();
	writeVTK ();
}
