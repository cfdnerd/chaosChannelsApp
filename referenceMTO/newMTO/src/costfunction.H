Info<<"cost function"<<endl;
MeanT=fvc::domainIntegrate(T).value()/area;
V=fvc::domainIntegrate(xh).value()/area-(voluse+0*Foam::max(0,0.5-opt*0.015));  

PowerDiss=0;
scalar PowerDissRef(0);
for (i=0; i<nObjPatch; i++)
{
     PowerDiss =PowerDiss - sum(
     phi.boundaryField()[conPatchList[i]]*(p.boundaryField()[conPatchList[i]]
     + 0.5*magSqr(U.boundaryField()[conPatchList[i]]))); 
}
PowerDissRef=PowerDiss;
PowerDiss=PowerDiss/PowerDiss0;
if(GeoDim==2)
{
 PowerDiss=PowerDiss/std::pow(mesh.V()[0],1.0/3);
 PowerDissRef=PowerDissRef/std::pow(mesh.V()[0],1.0/3);
}
reduce(PowerDiss, sumOp<scalar>()); 
reduce(PowerDissRef, sumOp<scalar>());

if (Pstream::master())
{ 
   ofstream outfile1("Voluse.txt",std::ios::app);
   outfile1<<V<<"\n";
   outfile1.close();

   ofstream outfile2("PowerDiss.txt",std::ios::app);
   if(opt==1)
   {
      outfile2<<"# Optimization Log: Power Dissipation Metrics\n";
      outfile2<<"# PowerDissNorm=constraint value (normalized); PowerDissRef=pre-normalization dissipation\n";
      outfile2<<"# Iteration PowerDissNorm PowerDissRef\n";
   }
   outfile2<<opt<<" "<<PowerDiss<<" "<<PowerDissRef<<"\n";
   outfile2.close(); 

   ofstream outfile3("MeanT.txt",std::ios::app);
   outfile3<<MeanT<<"\n";
   outfile3.close(); 

   ofstream outfile4("Time.txt",std::ios::app);
   outfile4<<runTime.elapsedClockTime()<<"\n";
   outfile4.close();    

   ofstream outfile5("maxT.txt",std::ios::app);
   outfile5<<maxT<<"\n";
   outfile5.close();
}
if(runTime.writeTime())
{
   T.write();
   //U.write();
}
