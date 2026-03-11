
 V=0;
 meanT=0;
 for (i=0;i<allcells;i++)
 {
     V=gamma[i]+V;
     meanT=meanT+T[i];
 }
 reduce(meanT, sumOp<scalar>());
 reduce(V, sumOp<scalar>());     
 V=V/nallcells-voluse;

 //DissPower
 scalar DissPowerRef(0);
 DissPower=0;
 for (i=0; i<nObjPatch; i++)
 {
      DissPower =DissPower - sum(
      phi.boundaryField()[conPatchList[i]]*(p.boundaryField()[conPatchList[i]]
      + 0.5*magSqr(U.boundaryField()[conPatchList[i]]))); 
 }
 DissPowerRef=DissPower;
 DissPower=DissPower/D_normalization;
 if(geo_dim==2) 
 {
   DissPower=DissPower/len;
   DissPowerRef=DissPowerRef/len;
 }
 reduce(DissPower, sumOp<scalar>()); 
 reduce(DissPowerRef, sumOp<scalar>());
 Info<<"constraint function (DissPower): "<<DissPower<<endl; 

 if (Pstream::master())
{
ofstream outfile1("Voluse.txt",std::ios::app);
outfile1<<V<<"\n";
outfile1.close();

ofstream outfile2("meanT.txt",std::ios::app);
outfile2<<meanT/nallcells<<"\n";
outfile2.close(); 
   
ofstream outfile3("PowerDiss.txt",std::ios::app);
if(opt==1)
{
outfile3<<"# Optimization Log: Power Dissipation Metrics\n";
outfile3<<"# PowerDissNorm=constraint value (normalized); PowerDissRef=pre-normalization dissipation\n";
outfile3<<"# Iteration PowerDissNorm PowerDissRef\n";
}
outfile3<<opt<<" "<<DissPower<<" "<<DissPowerRef<<"\n";
outfile3.close();

ofstream outfile4("Time.txt",std::ios::app);
outfile4<<runTime.elapsedClockTime()<<"\n";
outfile4.close(); 

}
