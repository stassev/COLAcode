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

////////////////////////////////////////////////////////////////////////////////////////////

/*
 * 
 * This file contains some standard functions for a PM code. Nothing COLA-specific.
 * 
 * 
 */ 




void WRt(float *d,int i,int j,int k,float f)
     {
     d[k + NGRID * (j + NGRID * i)]=f;
     }

float REd(float *d,int i,int j,int k)
     {
     return d[k + NGRID * (j + NGRID * i)];
     }



void WRtNROW(float *d,int i,int j,int k,float f)
     {
     d[k + NROW * (j + NROW * i)]=f;
     }



/////


  
void PtoMesh(float *densityz, struct particle_data *Pz)// Does CiC assignment.
      {
      
      int i;

     printf("Calculating PtoMesh\n");

     float scaleBox=((float) NGRID)/((float) BoxSize);
     float WPAR=pow(((float) NGRID)/((float)NROW),3);



     for(i=0; i<NGRID*NGRID*NGRID; i++)
            densityz[i] = -1.0;
            
     float X,Y,Z;
     float T1,T2,T3;
     float D1,D2,D3;
     float D2W,T2W;
     int iI,J,K,I1,J1,K1;

     for(i=1; i<=NumPart; i++)
     {
            X=Pz[i].Pos[0]*scaleBox;
            Y=Pz[i].Pos[1]*scaleBox;
            Z=Pz[i].Pos[2]*scaleBox;
            
            iI=(int) X;
            J=(int) Y;
            K=(int) Z;
            D1=X-((float) iI);
            D2=Y-((float) J);
            D3=Z-((float) K);
            T1=1.-D1;
            T2=1.-D2;
            T3=1.-D3;


            T2W =T2*WPAR;
            D2W =D2*WPAR;
            
//

            if(iI >= NGRID)iI=0;
            if(J >= NGRID)J=0;
            if(K >= NGRID)K=0;
            
            I1=iI+1;
               if(I1 >= NGRID)I1=0;
            J1=J+1;
               if(J1 >= NGRID)J1=0;
            K1=K+1;
               if(K1 >= NGRID)K1=0;


              WRt(densityz,iI ,J ,K ,REd(densityz,iI ,J ,K ) +T3*T1*T2W);
              WRt(densityz,I1,J ,K , REd(densityz,I1,J ,K ) +T3*D1*T2W);
              WRt(densityz,iI ,J1,K ,REd(densityz,iI ,J1,K ) +T3*T1*D2W);
              WRt(densityz,I1,J1,K , REd(densityz,I1,J1,K ) +T3*D1*D2W);
    
              WRt(densityz,iI ,J ,K1 ,REd(densityz,iI ,J ,K1 ) +D3*T1*T2W);
              WRt(densityz,I1,J ,K1 , REd(densityz,I1,J ,K1 ) +D3*D1*T2W);
              WRt(densityz,iI ,J1,K1 ,REd(densityz,iI ,J1,K1 ) +D3*T1*D2W);
              WRt(densityz,I1,J1,K1 , REd(densityz,I1,J1,K1 ) +D3*D1*D2W);
              
    }

      printf("CIC density assignment finished\n");

      return;
}


////////////////////////////////////////////////////////////////////////////////////////

void MtoParticles(struct particle_data *P, float *particles,float*mesh)// Does 3-linear interpolation
      {



     float X,Y,Z;
     float T1,T2,T3;
     float D1,D2,D3;
     float D2W,T2W;
     int iI,J,K,I1,J1,K1;
     int i;
     
     float scaleBox=((float) NGRID)/((float) BoxSize);
     float WPAR=1; 
     
     printf("Calculating MP...\n");
     
     
          for(i=1; i<=NumPart; i++)
     {
            X=P[i].Pos[0]*scaleBox;
            Y=P[i].Pos[1]*scaleBox;
            Z=P[i].Pos[2]*scaleBox;
            
            iI=(int) X;
            J=(int) Y;
            K=(int) Z;
            D1=X-((float) iI);
            D2=Y-((float) J);
            D3=Z-((float) K);
            T1=1.-D1;
            T2=1.-D2;
            T3=1.-D3;
            T2W =T2*WPAR;
            D2W =D2*WPAR;
//
            if(iI >= NGRID)iI=0;
            if(J >= NGRID)J=0;
            if(K >= NGRID)K=0;
            
            I1=iI+1;
               if(I1 >= NGRID)I1=0;
            J1=J+1;
               if(J1 >= NGRID)J1=0;
            K1=K+1;
               if(K1 >= NGRID)K1=0;



            particles[i-1] =                
              REd(mesh,iI,J,K    ) *T3*T1*T2W+
              REd(mesh,I1,J ,K   ) *T3*D1*T2W+
              REd(mesh,iI ,J1,K  ) *T3*T1*D2W+
              REd(mesh,I1,J1,K   ) *T3*D1*D2W+
              REd(mesh,iI ,J ,K1 ) *D3*T1*T2W+
              REd(mesh,I1,J ,K1  ) *D3*D1*T2W+
              REd(mesh,iI ,J1,K1 ) *D3*T1*D2W+
              REd(mesh,I1,J1,K1  ) *D3*D1*D2W;
              
          }
     
        
}

float putInBox(float din){
    float dout;
          dout=fmod(din,BoxSize);
          if (dout<0) dout+=BoxSize;   
    return dout;      
    }



void fft(float *arr,fftwf_complex *fft,int NROW)
{
    
fftwf_plan p=fftwf_plan_dft_r2c_3d(NROW,NROW,NROW,arr,fft,FFTW_ESTIMATE);

fftwf_execute(p);
fftwf_destroy_plan(p);
}
   




float meanM(float *arr,int N,float mean){ // Set mean of array values to mean.
        int i;
        float summ=0;
        float old;
        for(i=0;i<N;i++)summ+=arr[i];
        summ/=((float) N);
        old=summ;
        summ=mean-summ;
        for(i=0;i<N;i++)arr[i]+=summ;
        return old;
    }







void forces(fftwf_complex *P3D,int filter,int NGRID)
{

        N11=malloc(NGRID*NGRID*NGRID*sizeof(float));
        N12=malloc(NGRID*NGRID*NGRID*sizeof(float));
        N13=malloc(NGRID*NGRID*NGRID*sizeof(float));
        FN11=(fftwf_complex*) fftwf_malloc(NGRID*NGRID*((NGRID/2)+1) * sizeof(fftwf_complex));
        FN12=(fftwf_complex*) fftwf_malloc(NGRID*NGRID*((NGRID/2)+1) * sizeof(fftwf_complex));
        FN13=(fftwf_complex*) fftwf_malloc(NGRID*NGRID*((NGRID/2)+1) * sizeof(fftwf_complex));
        p11=fftwf_plan_dft_c2r_3d(NGRID,NGRID,NGRID,FN11,N11,FFTW_MEASURE);
        p12=fftwf_plan_dft_c2r_3d(NGRID,NGRID,NGRID,FN12,N12,FFTW_MEASURE);
        p13=fftwf_plan_dft_c2r_3d(NGRID,NGRID,NGRID,FN13,N13,FFTW_MEASURE);

 
  
    FN11[0]=0.0;
    FN12[0]=0.0;
    FN13[0]=0.0;
    

    
    int i,j;


    int iI, J,K,KMIN;
    float RK;
    complex float di,dj,dk,dens;



    float KK;


        for (iI=0;iI<NGRID/2+1;iI++){
        for (J=0;J<NGRID/2+1;J++){
            KMIN=0;
            if((iI==0)&&(J==0))KMIN=1;
        for (K=KMIN;K<(NGRID/2+1);K++){
                RK =(float) (K*K+iI*iI+J*J) ;
                
                di=(complex float) iI;
                dj=(complex float) J;
                dk=(complex float) K;
                dens = - P3D[K + (NGRID/2+1) * (J + NGRID * iI)];
                dens *= 1.0/RK/(pow((float) NGRID,3));
                KK=1;
                if (filter==1)KK=RK*Scale*Scale;
               
                    FN11[K + (NGRID/2+1) * (J + NGRID * iI)] = dens*I*di/Scale*KK;
                    FN12[K + (NGRID/2+1) * (J + NGRID * iI)] = dens*I*dj/Scale*KK;
                    FN13[K + (NGRID/2+1) * (J + NGRID * iI)] = dens*I*dk/Scale*KK;
                

            if ((J!=NGRID/2)&&(J!=0)){
                 j=NGRID-J;
                di=(complex float) iI;
                dj=-((complex float) J);
                dk=(complex float) K;
                dens = - P3D[K + (NGRID/2+1) * (j + NGRID * iI)];
                dens *= 1.0/RK/(pow((float) NGRID,3));
                    FN11[K + (NGRID/2+1) * (j + NGRID * iI)] = dens*I*di/Scale*KK;
                    FN12[K + (NGRID/2+1) * (j + NGRID * iI)] = dens*I*dj/Scale*KK;
                    FN13[K + (NGRID/2+1) * (j + NGRID * iI)] = dens*I*dk/Scale*KK;
                
            }
            
            if ((iI!=NGRID/2)&&(iI!=0)){
                 i=NGRID-iI;
                di=-((complex float) iI);
                dj=(complex float) J;
                dk=(complex float) K;
                dens = - P3D[K + (NGRID/2+1) * (J + NGRID * i)];

                dens *= 1.0/RK/(pow((float) NGRID,3));
            
                    FN11[K + (NGRID/2+1) * (J + NGRID * i)] = dens*I*di/Scale*KK;
                    FN12[K + (NGRID/2+1) * (J + NGRID * i)] = dens*I*dj/Scale*KK;
                    FN13[K + (NGRID/2+1) * (J + NGRID * i)] = dens*I*dk/Scale*KK;
                
            }
            
            if ((iI!=NGRID/2)&&(J!=NGRID/2)&&(J!=0)&&(iI!=0)){
                 i=NGRID-iI;
                 j=NGRID-J;
                di=-((complex float) iI);
                dj=-((complex float) J);
                dk=(complex float) K;
                dens = - P3D[K + (NGRID/2+1) * (j + NGRID * i)];

                dens *= 1.0/RK/(pow((float) NGRID,3));
                    FN11[K + (NGRID/2+1) * (j + NGRID * i)] = dens*I*di/Scale*KK;
                    FN12[K + (NGRID/2+1) * (j + NGRID * i)] = dens*I*dj/Scale*KK;
                    FN13[K + (NGRID/2+1) * (j + NGRID * i)] = dens*I*dk/Scale*KK;
                
            }
            
        }}}


      
    printf("fft ...\n");
      
      
         
         
         
    fftwf_execute(p11);
    fftwf_execute(p12);
    fftwf_execute(p13);
    printf("done\n");
// FINALIZE:


        fftwf_destroy_plan(p11);
        fftwf_destroy_plan(p12);
        fftwf_destroy_plan(p13);
        fftwf_free(FN11);
        fftwf_free(FN12);
        fftwf_free(FN13);
        
        // return;
    

// flush the output
    fflush(stdout) ;

}



void readSnapshots(void){


    struct particle_data *Ppointer[1];
    
    int    files;
    files=1;                               /* number of files per snapshot. CHANGE THIS AS NEEDED! */ 

    load_snapshot(input_fname,  files, Ppointer);
    P=Ppointer[0];
    
    NROW=(int) (0.5+powf(((float )NumPart),1./3.));
    NGRID=GridScale*NROW;
    
    printf("%i %i\n",NROW,NGRID);   

    Scale=2.*M_PI/BoxSize;

}






void ReconstructNormalOrder(float *dXr,float *dYr,float *dZr,struct particle_data * PR){
    //Starting with particle displacements arrays, calculate particle positions in a box.
    //Uses FORTRAN ordering
    
    int i,j,k,it;
    float dd;
    float ingrid=1./((float)NROW);
    

    
    for (i=0;i<NROW;i++){
        for (j=0;j<NROW;j++){
            for (k=0;k<NROW;k++){
              it=i+j*NROW+k*NROW*NROW+1;
                        dd=dXr[it-1];
          PR[it].Pos[0]=putInBox(((float)i)*ingrid*BoxSize+dd);
          
                        dd=dYr[it-1];
          PR[it].Pos[1]=putInBox(((float)j)*ingrid*BoxSize+dd);
          
                        dd=dZr[it-1];
          PR[it].Pos[2]=putInBox(((float)k)*ingrid*BoxSize+dd);
          
    }}}
    

}





void slice(struct particle_data *P){
    int it;
    printf("slice Begin\n");
    for(it=0; it<NumPart; it++){
      if (P[it+1].Pos[0]<BoxSize/32.){
          printf("%g %g %g\n",P[it+1].Pos[0],P[it+1].Pos[1],P[it+1].Pos[2]);
      }
    }
    printf("slice End\n");
    
}







void rearrange(float *dX,float*dY,float*dZ){ // convert between C ordering and fortran IC code ordering for gridded IC.
    int i,j,k,it,it1;
    float *dXtemp=malloc(NGRID*NGRID*NGRID*sizeof(float));
    float *dYtemp=malloc(NGRID*NGRID*NGRID*sizeof(float));
    float *dZtemp=malloc(NGRID*NGRID*NGRID*sizeof(float));
    
    for (i=0;i<NumPart;i++){
        dXtemp[i]=dX[i];
        dYtemp[i]=dY[i];
        dZtemp[i]=dZ[i];
    }
    
     for (i=0;i<NROW;i++){
          for (j=0;j<NROW;j++){
              for (k=0;k<NROW;k++){
                it=i+j*NROW+k*NROW*NROW;
                it1=k+j*NROW+i*NROW*NROW;

                dX[it]=dXtemp[it1];
                dY[it]=dYtemp[it1];
                dZ[it]=dZtemp[it1];

      }}}
    
    free(dXtemp);
    free(dYtemp);
    free(dZtemp);
}





////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
/*
 * 
 * The functions below are used when stepDistr == 2
 *   
 * 
 */
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////


double CosmoTimeFun (double a, void * params) {
       
       double f = a*a/Qfactor(a);
      
       return f;
     }
     
double   CosmoTime(double af)
     {
         // As defined below, 
         // CosmoTime(a)=int_0^a CosmoTimeFun(a) da
         // which corresponds to the usual time coordinate scaled by H0.
         
         double ai=1.0e-8;
       gsl_integration_workspace * w 
         = gsl_integration_workspace_alloc (5000);
       
       double result, error;
        double alpha=0;
     
       gsl_function F;
       F.function = &CosmoTimeFun;
       F.params = &alpha;
     
       gsl_integration_qag (&F, ai, af, 0, 1e-5, 5000,6,
                             w, &result, &error); 
     
       //printf ("result          = % .18f\n", result);
       //printf ("estimated error = % .18f\n", error);
    
       gsl_integration_workspace_free (w);
     
       return result;
     }

double AofTimeFun(double a,void * params){
    struct quadratic_params *p 
         = (struct quadratic_params *) params;
    float y = p->y;
    return CosmoTime(a)-y;
     }

float AofTime(float y) // Solves y=CosmoTime(a) for $a$
     {
       int status;
       int iter = 0, max_iter = 1000;
       const gsl_root_fsolver_type *T;
       gsl_root_fsolver *s;
       double r = 0;
       double x_lo = 0.0001, x_hi = 2.01;
       gsl_function F;
       //struct quadratic_params params = {1.0, 0.0, -5.0};
     
       struct quadratic_params params = {y};
     
       F.function = &AofTimeFun;
       F.params = &params;
     
       T = gsl_root_fsolver_brent;
       s = gsl_root_fsolver_alloc (T);
       gsl_root_fsolver_set (s, &F, x_lo, x_hi);
     
       //printf ("using %s method\n", 
               //gsl_root_fsolver_name (s));
     
       //printf ("%5s [%9s, %9s] %9s %9s\n",
               //"iter", "lower", "upper", "root", 
               //"err(est)");
     
       do
         {
           iter++;
           status = gsl_root_fsolver_iterate (s);
           r = gsl_root_fsolver_root (s);
           x_lo = gsl_root_fsolver_x_lower (s);
           x_hi = gsl_root_fsolver_x_upper (s);
           status = gsl_root_test_interval (x_lo, x_hi,
                                            0, 1.0e-6);
     
           //if (status == GSL_SUCCESS)
             //printf ("Converged:\n");
     
           //printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
                   //iter, x_lo, x_hi,
                   //r,
                   //x_hi - x_lo);
         }
       while (status == GSL_CONTINUE && iter < max_iter);
     
       gsl_root_fsolver_free (s);
     
       return r;
     }




