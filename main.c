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

    This is COLAcode: a serial particle mesh-based N-body code 
     illustrating the COLA (COmoving Lagrangian Acceleration) method 
     described in S. Tassev, M. Zaldarriaga, D. Eisenstein (2012).
     Check that paper (refered to as TZE below) for the details. 
     Before using the code make sure you read the README file as well as
     the Warnings section below.
    
    This version: Dec 18, 2012


*/

#include "stuff.h"


int main(int argc, char **argv)
{
    
    /*
    The code solves Newton's equation in the form (eq. (A.1) of TZE):
     d^2 (s)/dy^2 = - (3/2)*Omega_matter(a=1)*a*grad_x grad_x^(-2) \delta
    or its COLA modification. Above $y$ is defined as
    d/dy \equiv a/H_0*d/d\eta, 
    which is equivalent to T[]=d/dy in the notation of the paper.
    \eta = conformal time
    s = the displacement vector
    delta = fractional overdensity
    a = scale factor
    
    Velocity variable is defined as:
    vel\equiv T[s] = ds/dy =  Q*ds/da, 
    where Q\equiv a^3 H(a)/H0.
    
    When we subtract LPT, the new velocity variable is 
    d(s-s_LPT)/dy, so initially it is set to zero, as we see below.
    */
    
    
    /*
     Description of switches:
     
      useCOLA -- All it does is set sane values for the rest of the switches. 
         Note that "sane" values does not mean optimal!
         When set to 1 the code uses the COLA method. 
         Set to 0 for standard non-COLA N-body.
         
      subtractLPT -- Whether to subtract LPT displacements. 
         Set this to 1 to subtract the LPT displacements as the COLA method prescribes. 
         Otherwise, set to 0 for standard non-COLA approach.
         
      stepDistr -- How to distribute the timesteps between 
         the initial and final time.
         stepDistr=0 -- uniformly in $a$ 
         stepDistr=1 -- uniformly in $log(a)$. 
         stepDistr=2 -- uniformly distributed in a user defined
            variable E=CosmoTime(a). As an example E is set to cosmic time in the code.
      
      StdDA -- Whether to assume standard or non-standard time-dependence for 
         the (residual) displacement in the non-COLA (COLA) implementation.
         StdDA=0 is the modified COLA time-stepping (see eq. (A.15) in TZE).
         StdDA=1 is non-integral standard stepping. (Avoid! Just for debugging.)
         StdDA=2 is the integral (Quinn et al (1997), astro-ph/9710043) 
           standard stepping. Equivalent to equations A.3 (non-COLA) and A.7 (COLA) in TZE.
         
      When StdDA=0, one needs to set fullT and nLPT.
         fullT=0 assumes time dependence for velocity = A + B a^nLPT, with A>>B a^nLPT. (A and B are irrelevant)
         fullT=1 assumes time dep. for velocity = B a^nLPT
         nLPT is a real number. Sane values lie in the range (-4,3.5). Cannot be 0, but of course can be -> 0 (say 0.001).
         See Section A.3 of TZE.
      
     */
     
     
     
     /* !!!!!!!!!!!!!!!!!!!!!!!!!!
                WARNINGS 
      * !!!!!!!!!!!!!!!!!!!!!!!!!!
      
        0. We do not claim the code to be optimized in any way! 
           It is written more as an illustration of how to 
           implement the COLA method in an existing N-body code 
           rather than for production purposes.
       
        1. As written, the code assumes grid-like initial conditions. The reason
           lies in the way we extract the ZA and 2LPT displacements from the 2LPTic code output.
           This is easily fixable for glass-like IC but is not done below.
        
        2. The code assumes input from Roman Scoccimarro's serial FORTRAN 2LPTic code.
           Since that code is in FORTRAN, indices of flat 3D arrays are reversed ordered compared to C. 
           We use rearrange() function below to fix that. For other IC, one may need to 
           remove the calls to that function, as well as change the FORTRAN-style orderings in some of
           the other functions.
        
        3. Filenames, number of files per snapshot, OmegaMatter, z_{initial}, z_{final} 
           and switch values are all hardcoded below. Change them appropriately!
         
        4. As written, the input needded by the code is a ZA snapshot from 2LPTic at z=0, 
           and a second snapshot from 2LPTic for ZA+2LPT again at z=0. The initial 
           conditions at z_initial are calculated in the code below from those two snapshots.
            
        5. For the output of the code, check the comments in the end of main().   
        
        6. Assume flat cosmology throughout.
        
     */
     
     
     
     
     
    int stepDistr,nsteps;
    
    float zi,zf;        

    Om=0.2768; // Omega_CDM+Omega_baryon today.
    zi=9.; // Initial redshift. Works well for COLA.
    zf=0.0; // Final redshift. CHANGE THIS AS NEEDED!
    nsteps=10; // Number of time steps.
    GridScale=3; // How many times finer (in 1-dim) is the PM grid compared to the particle grid.

    printf("OmegaMatter(a=1) = %g\n",Om) ;
    printf("zi = %g\n",zi) ;
    printf("zf = %g\n",zf) ;
    printf("nsteps = %i\n",nsteps) ;
    printf("GridScale = %i\n",GridScale) ;





    
/* **************************
 * **************************
 * 
 * Set switch options.
 * 
 * **************************
 * **************************
 */ 
    
    
    
    
    int useCOLA=1; // Set "sane" options for using COLA or non-COLA

    
    if (useCOLA==1){
        subtractLPT=1; 
        stepDistr=0;
        StdDA=0; // Set to 0 if modified time-stepping is to be used. Set to 2 otherwise.
    }
    else{
        subtractLPT=0; 
        stepDistr=1;
        StdDA=2;
    }
    
    printf("subtractLPT = %g\n",subtractLPT) ;
    printf("stepDistr = %i\n",stepDistr) ;
    printf("StdDA = %i\n",StdDA) ;

    
    if (StdDA==0){
        fullT=1;
        nLPT=-2.5;
        printf("fullT = %i\n",fullT) ;
        printf("nLPT = %g\n",nLPT) ;
    }





/* **************************
* **************************
* 
* All switch options set.
* 
* **************************
* **************************
*/ 
    
    
    
    
    
    float Om143;
    float ai,af, A,da;
    // Initial and final scale factors:
    ai=1.0/(1.0+zi);
    af=1.0/(1.0+zf);
    
    A=ai; // This is the $a$ which we'll be advancing below.
    
    if (stepDistr==0) da=(af-ai)/((float)nsteps);
    if (stepDistr==1) da=(log(af)-log(ai))/((float)nsteps);
    if (stepDistr==2) da=(CosmoTime(af)-CosmoTime(ai))/((float)nsteps);





/* **************************
 * **************************
 * 
 * Import ZA and 2LPT displacements.
 * 
 * **************************
 * **************************
 */ 





    // Import ZA+2LPT output from 2LPTic code at z=0.
    sprintf(input_fname,"/scratch/stassev/G_FULL/G0_100Mpc/proba57lpt2");    // CHANGE!
    readSnapshots(); 
        // Parameters imported from snapshot:
        // NROW = (NumPart)^(1/3) is 1-dim particle grid size.
        // NGRID=GridScale*NROW is PM grid size
        // BoxSize in [Mpc/h]
    Displacements(); // Extract displacements from snapshot.

    float *dX2=malloc(NumPart*sizeof(float));
    float *dY2=malloc(NumPart*sizeof(float));
    float *dZ2=malloc(NumPart*sizeof(float));

    int i;
    for(i=0; i<NumPart; i++)
    {
        dX2[i]=dX[i];
        dY2[i]=dY[i];
        dZ2[i]=dZ[i];
    }
    rearrange(dX2,dY2,dZ2); // rearrange indices from FORTRAN -> C!
    free(dX);
    free(dY);
    free(dZ);
        

    // Import ZA output from 2LPTic code at z=0.
    sprintf(input_fname,"/scratch/stassev/G_FULL/G0_100Mpc/proba57za");    // CHANGE!
    readSnapshots();  
    Displacements();  // Extract displacements from snapshot.
    rearrange(dX,dY,dZ);
    
    float *dXz=malloc(NumPart*sizeof(float));
    float *dYz=malloc(NumPart*sizeof(float));
    float *dZz=malloc(NumPart*sizeof(float));

    float Di=growthD(A); // initial growth factor
    float Di2=growthD2(A); // initial 2nd order growth factor
    

    for(i=0; i<NumPart; i++)
    {
        dX2[i]-=dX[i]; // this stores 2nd-order LPT displ.
        dY2[i]-=dY[i];
        dZ2[i]-=dZ[i];
        
        dXz[i]=dX[i]; // this stores ZA displ.
        dYz[i]=dY[i];
        dZz[i]=dZ[i];
    }

    
    
    float Dv=DprimeQ(A,1.0); // T[D_{za}]=dD_{za}/dy
    float Dv2=growthD2v(A);  // T[D_{2lpt}]=dD_{2lpt}/dy
    
    for(i=0; i<NumPart; i++)
    {        
        dX[i]=dXz[i]*Di+dX2[i]*Di2;
        dY[i]=dYz[i]*Di+dY2[i]*Di2;
        dZ[i]=dZz[i]*Di+dZ2[i]*Di2;
        
        if (subtractLPT <0.5){ // If subtractLPT=0 (non-COLA), then velocity is ds/dy, which is simply the 2LPT IC.
            P[i+1].Vel[0]=dXz[i]*Dv+dX2[i]*Dv2;
            P[i+1].Vel[1]=dYz[i]*Dv+dY2[i]*Dv2;
            P[i+1].Vel[2]=dZz[i]*Dv+dZ2[i]*Dv2;
        }
        else { // Set vel=0 if we subtract LPT. 
               //This is the same as the action of the operator L_- from TZE, as initial velocities are in 2LPT.
            P[i+1].Vel[0]=0.0;
            P[i+1].Vel[1]=0.0;
            P[i+1].Vel[2]=0.0;
        }
    }

    // Set P[i].Pos[j] array in Gadget format from displacements arrays.
    ReconstructNormalOrder(dX,dY,dZ,P); 
    
    // Now P holds the particle structure containing the initial 
    // particle positions and velocities as obtained in 2LPT. 
    // Note that the initial velocities are set to 0 in the COLA method,
    // due to the action of the L_- operator in TZE.
    
    
    

/* **************************
 * **************************
 * 
 * Importing ZA and 2LPT displacements is done.
 * 
 * **************************
 * **************************
 */ 


    
    density=malloc(NGRID*NGRID*NGRID*sizeof(float)); // PM density grid.
    P3D = (fftwf_complex*) fftwf_malloc(NGRID*NGRID*((NGRID/2)+1) * sizeof(fftwf_complex));  //will contain fftw'ed density
    
    float growth1=Di;
    float growth1L2=Di2;
    
    float dyyy,dda;
    float ax,ay,az;
    
    float q1,q2;
    float da2;

    float da1;
    float AF,AI,AC,AFF;
    float sumx,sumy,sumz;
    
    printf("Init done \n");
    
    
    
    
/* **************************
 * **************************
 * 
 * Evolve with time using KDK!
 * 
 * **************************
 * **************************
 */ 
 
 
 
    
    int timeStep;
    
    // AI stores the scale factor to which the velocities have been kicked to. Initially it's just A.
    AI=A;
for (timeStep=0;timeStep<nsteps+1;timeStep++){
    
    //AFF is the scale factor to which we should drift the particle positions.
    //AF is the scale factor to which we should kick the particle velocities.
    
    if (stepDistr==0) AFF=A+da;
    if (stepDistr==1) AFF=A*exp(da);
    if (stepDistr==2) AFF=AofTime(CosmoTime(A)+da);

    if (timeStep==nsteps) AF=A;// half time-step for final kick
    else { // Set to mid-point of interval. How one chooses the mid-point is really an extra degree of freedom in the code. But I find negligible effects from the different choices below. So, this is not exported as an extra switch at this point.
        if (stepDistr==0) AF=A+da*0.5;
        if (stepDistr==1) AF=A*exp(da*0.5);
        if (stepDistr==2) AF=AofTime((CosmoTime(AFF)+CosmoTime(A))*0.5); 
    }
    sumx=0; sumy=0; sumz=0;
    printf("iteration = %i begin \n",timeStep);    fflush(stdout);
    
/*
 * 
 *  Force calculation
 * 
 */ 
    
    PtoMesh(density, P);// Do Cloud-in-Cell assignment. Maps P to density.
    fft(density,P3D,NGRID);//P3D=fftw'ed density
    
    //Next calculate acceleration from fftw'ed density. 
    //This returns N11,N12,N13 which hold the components of
    // the vector (grad grad^{-2} density) on a grid.
    forces(P3D,0,NGRID); 
    
    //Now find the accelerations at the particle positions using 3-linear interpolation. 
    // Result is stored in the (dX,dY,dZ) arrays.
    MtoParticles(P, dX,N11); 
    MtoParticles(P, dY,N12);
    MtoParticles(P, dZ,N13);
    free(N11);
    free(N12);
    free(N13);
    
    //Make sure we'll conserve momentum.
    meanM(dX,NumPart,0.);
    meanM(dY,NumPart,0.);
    meanM(dZ,NumPart,0.);


/*
 * 
 *  KICK:
 * 
 */ 

    Om143=pow(Om/(Om+(1-Om)*A*A*A),1./143.);
    
    if      (StdDA==0)  dda=Sphi(AI,AF,A);
    else if (StdDA==1)  dda=(AF-AI)*A/Qfactor(A);
    else                dda=SphiStd(AI,AF);
    

    
    q2=1.5*Om*growth1*growth1*(1.0+7./3.*Om143)*A; // T^2[D_{2lpt}]=d^2 D_{2lpt}/dy^2
    q1=1.5*Om*growth1*A; // T^2[D_{ZA}]=d^2 D_{ZA}/dy^2
    

    for(i=0; i<NumPart; i++)
     {
        ax=-1.5*Om*dX[i]-subtractLPT*(dXz[i]*q1+dX2[i]*q2)/A;
        ay=-1.5*Om*dY[i]-subtractLPT*(dYz[i]*q1+dY2[i]*q2)/A;
        az=-1.5*Om*dZ[i]-subtractLPT*(dZz[i]*q1+dZ2[i]*q2)/A;
        
        P[i+1].Vel[0] += ax*dda;
        P[i+1].Vel[1] += ay*dda;
        P[i+1].Vel[2] += az*dda;

        sumx+=P[i+1].Vel[0];
        sumy+=P[i+1].Vel[1];
        sumz+=P[i+1].Vel[2];
    }
    

    
    fflush(stdout);
    
    sumx/=(float)NumPart; // will subtract this below to conserve momentum. Should be conserved ... but just in case 3-linear interpolation makes a problem. Never checked whether this makes a difference.
    sumy/=(float)NumPart;
    sumz/=(float)NumPart;
    
    if (timeStep==nsteps){
        //At final timestep, add back LPT velocities if we had subtracted them. This corresponds to L_+ operator in TZE.
        Dv=DprimeQ(A,1.0); // dD_{za}/dy
        Dv2=growthD2v(A); // dD_{2lpt}/dy
        for(i=0; i<NumPart; i++)
        {
            P[i+1].Vel[0]+=-sumx+(dXz[i]*Dv+dX2[i]*Dv2)*subtractLPT;
            P[i+1].Vel[1]+=-sumy+(dYz[i]*Dv+dY2[i]*Dv2)*subtractLPT;
            P[i+1].Vel[2]+=-sumz+(dZz[i]*Dv+dZ2[i]*Dv2)*subtractLPT;
        }
        goto finalize; // Sorry for "goto" :)
    }
    
    
/*
 * 
 * DRIFT:
 * 
 */ 
    
    //keep in mind this renaming:
    AC=AF;
    AF=AFF;
    
    if (StdDA == 0)    dyyy=Sq(A,AF,AC);
    else if (StdDA==1) dyyy=(AF-A)/Qfactor(AC);
    else               dyyy=SqStd(A,AF);

    
    da1=growthD(AF)-growth1; // change in D
    da2=growthD2(AF)-growth1L2; // change in D_{2lpt}
    

    for(i=0; i<NumPart; i++)
     {
        P[i+1].Pos[0] +=(P[i+1].Vel[0]-sumx)*dyyy;
        P[i+1].Pos[1] +=(P[i+1].Vel[1]-sumy)*dyyy;
        P[i+1].Pos[2] +=(P[i+1].Vel[2]-sumz)*dyyy;

        P[i+1].Pos[0]=putInBox(P[i+1].Pos[0]+subtractLPT*(dXz[i]*da1+dX2[i]*da2));
        P[i+1].Pos[1]=putInBox(P[i+1].Pos[1]+subtractLPT*(dYz[i]*da1+dY2[i]*da2));
        P[i+1].Pos[2]=putInBox(P[i+1].Pos[2]+subtractLPT*(dZz[i]*da1+dZ2[i]*da2));    
    }
    
/*
 * 
 * STEP IN TIME:
 * 
 */
  
    A=AF;
    AI=AC;

  
    growth1=growthD(A);
    growth1L2=growthD2(A);
    
    printf("iteration = %i done \n",timeStep);
    printf("a = %g \n",A);
    printf("z = %g \n",1.0/A-1.0);
    fflush(stdout);
}
    finalize:
    
    // At this point P contains the particle positions and velocities at redshift zf in Gadget units.
    
    printf("slice:\n"); // Output a slice just for the sake of doing something with P.
    slice(P);
    fflush(stdout);

   
   
// Surely I've missed to free a lot of the malloc'ed arrays ... 

    
    
    return 0;
}










///////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////////////////////////////////////////////////////////////////////////  
///////////////////////////////////////////////////////////////////////////////////////////  










/*
 * 
 * Extracting displacements from the 2LPTic snapshots.
 * 
 * 
 * 
 */ 
 

void Displacements(void){
    
    dX=malloc(NumPart*sizeof(float));
    dY=malloc(NumPart*sizeof(float));
    dZ=malloc(NumPart*sizeof(float));
   
    obtainDisplacements(P,dX,dY,dZ);
    
}
 

float fix(float x,float q){ // Extract particle displacement from Eulerian and Lagrangian position. It takes care of periodic BC assuming box size is >> rms particle displacement (~15Mpc at z=0).
    q=putInBox(q);
    float s=x-q;
    
    if (-s>BoxSize/2.)s=x-(q-BoxSize);
    if (s>BoxSize/2.) s=(x-BoxSize)-q;
    
    return s;
    
}

void obtainDisplacements(struct particle_data *P,float *dX,float *dY,float *dZ){ // uses FORTRAN ordering

      int i,j,k,it;
      float ingrid=1./((float)NROW);
    
      float summ,f;
      summ=0;

      for (i=0;i<NROW;i++){
          for (j=0;j<NROW;j++){
              for (k=0;k<NROW;k++){
                it=i+j*NROW+k*NROW*NROW+1;

                f=fix(P[it].Pos[0],((float)i)*ingrid*BoxSize);
                WRtNROW(dX, i, j, k, f);
                
                summ+=f;
                f=fix(P[it].Pos[1],((float)j)*ingrid*BoxSize);
                WRtNROW(dY, i, j, k, f);

                f=fix(P[it].Pos[2],((float)k)*ingrid*BoxSize);
                WRtNROW(dZ, i, j, k, f);

      }}}
    // should be ->zero:
    //printf("==============AVERAGE 1D DISPL = %f\n",summ*powf(ingrid,3));
    meanM(dX,NROW*NROW*NROW,0.0);
    meanM(dY,NROW*NROW*NROW,0.0);
    meanM(dZ,NROW*NROW*NROW,0.0);
    
    
  }




/*
 * 
 *  Below follow the time dependent functions used in the code.
 *  Start with growth factors and derivatives.
 * 
 */



float growthD(float a){ // growth factor for LCDM
    return growthDtemp(a)/growthDtemp(1.0);
}


float growthDtemp(float a){
    // Decided to use the analytic expression for LCDM. More transparent if I change this to numerical integration?
    float x=-Om/(Om-1.0)/(a*a*a);

    
    
    float hyperP=0,hyperM=0;
    
    if (fabs(x-1.0) < 1.e-3) {
       // printf("mechka\n");
    hyperP=0.859596768064608 - 0.1016599912520404*(-1.0 + x) + 0.025791094277821357*pow(-1.0 + x,2) - 0.008194025861121475*pow(-1.0 + x,3) + 0.0029076305993447644*pow(-1.0 + x,4) - 0.0011025426387159761*pow(-1.0 + x,5) + 0.00043707304964624546*pow(-1.0 + x,6) - 0.0001788889964687831*pow(-1.0 + x,7);
    hyperM=1.1765206505266006 + 0.15846194123099624*(-1.0 + x) - 0.014200487494738975*pow(-1.0 + x,2) + 0.002801728034399257*pow(-1.0 + x,3) - 0.0007268267888593511*pow(-1.0 + x,4) + 0.00021801569226706922*pow(-1.0 + x,5) - 0.00007163321597397065*pow(-1.0 + x,6) +    0.000025063737576245116*pow(-1.0 + x,7);
    }
    else {
        if (x < 1.0) {
            hyperP=gsl_sf_hyperg_2F1(1.0/2.0,2.0/3.0,5.0/3.0,-x);
            hyperM=gsl_sf_hyperg_2F1(-1.0/2.0,2.0/3.0,5.0/3.0,-x);
        }
        x=1.0/x;
        if ((x < 1.0) && (x>1.0/30)) {
            
            hyperP=gsl_sf_hyperg_2F1(-1.0/6.0,0.5,5.0/6.0,-x);
            hyperP*=4*sqrt(x);
            hyperP+=-3.4494794123063873799*pow(x,2.0/3.0);
        
            hyperM=gsl_sf_hyperg_2F1(-7.0/6.0,-0.5,-1.0/6.0,-x);
            hyperM*=4.0/7.0/sqrt(x);
            hyperM+=pow(x,2.0/3.0)*(-1.4783483195598803057); //-(Gamma[-7/6]*Gamma[5/3])/(2*sqrt[Pi])
        }
        if (x<=1.0/30.0){
            hyperP=3.9999999999999996*sqrt(x) - 3.4494794123063865*pow(x,0.6666666666666666) + 0.3999999999999999*pow(x,1.5) -    0.13636363636363635*pow(x,2.5) + 0.07352941176470587*pow(x,3.5) - 0.04755434782608695*pow(x,4.5) +    0.033943965517241374*pow(x,5.5) - 0.02578125*pow(x,6.5) + 0.020436356707317072*pow(x,7.5) -    0.01671324384973404*pow(x,8.5) + 0.013997779702240564*pow(x,9.5) - 0.011945562847590041*pow(x,10.5) + 0.01035003662109375*pow(x,11.5) - 0.009080577904069926*pow(x,12.5);
            hyperM=0.5714285714285715/sqrt(x) + 2.000000000000001*sqrt(x) - 1.4783483195598794*pow(x,0.66666666666666666) +    0.10000000000000002*pow(x,1.5) - 0.022727272727272735*pow(x,2.5) + 0.009191176470588237*pow(x,3.5) -    0.004755434782608697*pow(x,4.5) + 0.002828663793103449*pow(x,5.5) - 0.0018415178571428578*pow(x,6.5) +    0.0012772722942073172*pow(x,7.5) - 0.0009285135472074472*pow(x,8.5) + 0.0006998889851120285*pow(x,9.5) -    0.0005429801294359111*pow(x,10.5) + 0.0004312515258789064*pow(x,11.5) - 0.00034925299631038194*pow(x,12.5);
        }
    }
    
    
    if (a>0.2) return sqrt(1.0 + (-1.0 + pow(a,-3))*Om)*(3.4494794123063873799*pow(-1.0 + 1.0/Om,0.666666666666666666666666666) + (hyperP*(4*pow(a,3)*(-1.0 + Om) - Om) - 7.0*pow(a,3)*hyperM*(-1.0 + Om))/(pow(a,5)*(-1.0+ Om) - pow(a,2)*Om));
    return (a*pow(1 - Om,1.5)*(1291467969*pow(a,12)*pow(-1 + Om,4) + 1956769650*pow(a,9)*pow(-1 + Om,3)*Om + 8000000000*pow(a,3)*(-1 + Om)*pow(Om,3) + 37490640625*pow(Om,4)))/(1.5625e10*pow(Om,5));
    
        }






float Qfactor(float a){ // Q\equiv a^3 H(a)/H0.
    return sqrt(Om/(a*a*a)+1.0-Om)*a*a*a;
}




float growthD2(float a){// Second order growth factor
    return growthD2temp(a)/growthD2temp(1.0);
}


float growthD2temp(float a){
    float d= growthD(a);
    float omega=Om/(Om+(1.0-Om)*a*a*a);
    return d*d*pow(omega,-1./143.);
}
 
float growthD2v(float a){ // explanation is in main()
    float d2= growthD2(a);
    float omega=Om/(Om+(1.0-Om)*a*a*a);
    return Qfactor(a)*(d2/a)*2.0*pow(omega,6./11.);
}

 

float decayD(float a){ // D_{-}, the decaying mode
    return sqrt(Om/(a*a*a)+1.0-Om);
}



double DprimeQ(double a,float nGrowth){ // returns Q*d(D_{+}^nGrowth*D_{-}^nDecay)/da, where Q=Qfactor(a)
        float nDecay=0.0;// not interested in decay modes in this code.
        float Nn=6.0*pow(1.0 - Om,1.5)/growthDtemp(1.0);
        return (pow(decayD(a),-1.0 + nDecay)*pow(growthD(a),-1.0 + nGrowth)*(nGrowth*Nn- (3.0*(nDecay + nGrowth)*Om*growthD(a))/(2.*a)));
        
    }
   


/*
 * 
 *  
 * Functions for our modified time-stepping (used when StdDA=0):
 * 
 * 
 */


double fun (double a, void * params) {
       
       double f;
       if (fullT==1) f = gpQ(a)/Qfactor(a); 
       else f = 1.0/Qfactor(a);
      
       return f;
     }
     
double   Sq(double ai,double af,double aRef)
     {
       gsl_integration_workspace * w 
         = gsl_integration_workspace_alloc (5000);
       
       double result, error;
        double alpha=0;
     
       gsl_function F;
       F.function = &fun;
       F.params = &alpha;
     
       gsl_integration_qag (&F, ai, af, 0, 1e-5, 5000,6,
                             w, &result, &error); 
     
       //printf ("result          = % .18f\n", result);
       //printf ("estimated error = % .18f\n", error);
      
     
       gsl_integration_workspace_free (w);
     
      
        //   result  = \int fun(a)da 
      
       if (fullT==1) return result/gpQ(aRef);
       else return result;
       
     }
     
     
double   Sphi(double ai,double af,double aRef)
     {
       double result;
       result=(gpQ(af)-gpQ(ai))*aRef/Qfactor(aRef)/DERgpQ(aRef);
     
     
       return result;
     }



double gpQ(double a){ 
    
        return pow(a,nLPT);
        
    }
double DERgpQ(double a){ // This must return d(gpQ)/da
        
        return nLPT*pow(a,nLPT-1);
        
    }
    




/*
 * 
 *  
 * Functions for Quinn et al time-stepping (used when StdDA=2):
 * 
 * 
 */ 


double funSqStd (double a, void * params) {
       
       double f = 1.0/Qfactor(a);
      
       return f;
     }
     
double   SqStd(double ai,double af)
     {
       gsl_integration_workspace * w 
         = gsl_integration_workspace_alloc (5000);
       
       double result, error;
        double alpha=0;
     
       gsl_function F;
       F.function = &funSqStd;
       F.params = &alpha;
     
       gsl_integration_qag (&F, ai, af, 0, 1e-5, 5000,6,
                             w, &result, &error); 
     
       //printf ("result          = % .18f\n", result);
       //printf ("estimated error = % .18f\n", error);
      
     
       gsl_integration_workspace_free (w);
     
      
       return result;
       
       
     }
     
     

double funSphiStd (double a, void * params) {
       
       double f = a/Qfactor(a);
      
       return f;
     }
     
double   SphiStd(double ai,double af)
     {
       gsl_integration_workspace * w 
         = gsl_integration_workspace_alloc (5000);
       
       double result, error;
        double alpha=0;
     
       gsl_function F;
       F.function = &funSphiStd;
       F.params = &alpha;
     
       gsl_integration_qag (&F, ai, af, 0, 1e-5, 5000,6,
                             w, &result, &error); 
     
       //printf ("result          = % .18f\n", result);
       //printf ("estimated error = % .18f\n", error);
      
     
       gsl_integration_workspace_free (w);
     
       return result;
       
       
     }





