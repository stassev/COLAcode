#include "stuff.h"

/*
    Copyright (c) 2011-2013       Svetlin Tassev
                           Harvard University, Princeton University
 
    This file is part of COLAcode.

    COLAcode is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    COLAcode is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with COLAcode.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 * Copyright (c) 2005       Volker Springel
 *                          Max-Plank-Institute for Astrophysics
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * 
 * 
 *  INSERTING SOME GADGET-2 CODE TO READ SNAPSHOTS:
 *  (slightly modified to set pointers and global variables)
 * 
 * 
 */ 


/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int load_snapshot(char *fname, int files,  struct particle_data **Ppointer)
{
  FILE *fd;
  char   buf[200];
  int    i,k,dummy,ntot_withmasses;
  int    n,pc,pc_new,pc_sph;
struct io_header_1 header1;
struct particle_data *P;
int *Id;
int *IdPointer[1];

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i=0, pc=1; i<files; i++, pc=pc_new)
    {
      if(files>1)
    sprintf(buf,"%s.%d",fname,i);
      else
    sprintf(buf,"%s",fname);

      if(!(fd=fopen(buf,"r")))
    {
      printf("can't open file `%s`\n",buf);
      exit(0);
    }

      printf("reading `%s' ...\n",buf); fflush(stdout);

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);



      if(files==1)
    {
      for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
        NumPart+= header1.npart[k];
    //  Ngas= header1.npart[0];
    }
      else
    {
      for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
        NumPart+= header1.npartTotal[k];
    //  Ngas= header1.npartTotal[0];
    }

      for(k=0, ntot_withmasses=0; k<5; k++)
    {
      if(header1.mass[k]==0)
        ntot_withmasses+= header1.npart[k];
    }


/////////////////////////
    if(i==0){
    allocate_memory(IdPointer,Ppointer);
    Id=IdPointer[0];
    P=Ppointer[0];
      }
/////////////////////////

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
    {
      for(n=0;n<header1.npart[k];n++)
        {
          fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
          pc_new++;
        }
    }
      SKIP;

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
    {
      for(n=0;n<header1.npart[k];n++)
        {
          fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
          pc_new++;
        }
    }
      SKIP;
    

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
    {
      for(n=0;n<header1.npart[k];n++)
        {
          fread(&Id[pc_new], sizeof(int), 1, fd);
          pc_new++;
        }
    }
      SKIP;


      if(ntot_withmasses>0)
    SKIP;
      for(k=0, pc_new=pc; k<6; k++)
    {
      for(n=0;n<header1.npart[k];n++)
        {
          P[pc_new].Type=k;

          if(header1.mass[k]==0)
        fread(&P[pc_new].Mass, sizeof(float), 1, fd);
          else
        P[pc_new].Mass= header1.mass[k];
          pc_new++;
        }
    }
      if(ntot_withmasses>0)
    SKIP;
      

      if(header1.npart[0]>0)
    {
      SKIP;
      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
        {
          fread(&P[pc_sph].U, sizeof(float), 1, fd);
          pc_sph++;
        }
      SKIP;

      SKIP;
      for(n=0, pc_sph=pc; n<header1.npart[0];n++)
        {
          fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
          pc_sph++;
        }
      SKIP;

      if(header1.flag_cooling)
        {
          SKIP;
          for(n=0, pc_sph=pc; n<header1.npart[0];n++)
        {
          fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
          pc_sph++;
        }
          SKIP;
        }
      else
        for(n=0, pc_sph=pc; n<header1.npart[0];n++)
          {
        P[pc_sph].Ne= 1.0;
        pc_sph++;
          }
    }

      fclose(fd);
    }


  reordering(Id, P);  /* call this routine only if your ID's are set properly */
  BoxSize=header1.BoxSize;
  printf("Redshift = %g\n",header1.redshift);

return 0;
}

///////////////////////////////////////////////


/* this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(int**IdPointer, struct particle_data **Ppointer)
{
      struct particle_data *P;
      int * Id;
  printf("allocating memory...\n");

  if(!(P=malloc(NumPart*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  P--;   /* start with offset 1 */

  
  if(!(Id=malloc(NumPart*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  Id--;   /* start with offset 1 */




  Ppointer[0]=P;
  IdPointer[0]=Id;

  printf("allocating memory...done\n");
  return 0;
}

///////////////////////////////////////////////


/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 */
int reordering(int*Id, struct particle_data *P)
{
  int i;
  int idsource, idsave, dest;
  struct particle_data psave, psource;


  printf("reordering....\n");

  for(i=1; i<=NumPart; i++)
    {
      if(Id[i] != i)
    {
      psource= P[i];
      idsource=Id[i];
      dest=Id[i];

      do
        {
          psave= P[dest];
          idsave=Id[dest];

          P[dest]= psource;
          Id[dest]= idsource;
          
          if(dest == i) 
        break;

          psource= psave;
          idsource=idsave;

          dest=idsource;
        }
      while(1);
    }
    }

  printf("done.\n");

  Id++;   
  free(Id);

  printf("space for particle ID freed\n");
return 0;
}

//
// END OF GADGET BLOB.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
