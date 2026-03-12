
scalar diff(volScalarField &gamma,const scalarField &V,double del,double eta,int n)
{
     int i;
     scalar z=0;
     double *x =new double[n];
     const scalar etaEps(1e-6);
     const scalar etaSafe = Foam::min(Foam::max(eta, etaEps), scalar(1.0)-etaEps);
     const scalar etaDenominator = Foam::max(etaSafe, etaEps);
     const scalar oneMinusEtaDenominator = Foam::max(scalar(1.0)-etaSafe, etaEps);
     
     for(i=0;i<n;i++)
     {
        const scalar gammaClamped = Foam::min(Foam::max(gamma[i], scalar(0.0)), scalar(1.0));
        if(gammaClamped<=etaSafe)
        {
          const scalar ratio = gammaClamped/etaDenominator;
          x[i]=etaSafe*(Foam::exp(-del*(1-ratio))-(1-ratio)*Foam::exp(-del));
        }
        else
        {
          const scalar ratio = (gammaClamped-etaSafe)/oneMinusEtaDenominator;
          x[i]=etaSafe+(1-etaSafe)*(1-Foam::exp(-del*ratio)+ratio*Foam::exp(-del));
        }
     }
     for(i=0;i<n;i++)
     {
        z=z+(gamma[i]-x[i])*V[i];
     }
     delete[] x;
     return {z};
}
