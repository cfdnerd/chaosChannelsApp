
scalar diff
(
  volScalarField &gamma,
  const scalarField &V,
  double del,
  double eta,
  int n
)
{
     int i;
     scalar z=0;
     double *x =new double[n];
     const scalar etaProjection(eta);
     const scalar etaDenominator(etaProjection);
     const scalar oneMinusEtaDenominator(scalar(1.0)-etaProjection);
     
     for(i=0;i<n;i++)
     {
        const scalar gammaProjection = gamma[i];
        if(gammaProjection<=etaProjection)
        {
          const scalar ratio = gammaProjection/etaDenominator;
          x[i]=etaProjection*(Foam::exp(-del*(1-ratio))-(1-ratio)*Foam::exp(-del));
        }
        else
        {
          const scalar ratio = (gammaProjection-etaProjection)/oneMinusEtaDenominator;
          x[i]=etaProjection+(1-etaProjection)*(1-Foam::exp(-del*ratio)+ratio*Foam::exp(-del));
        }
     }
     for(i=0;i<n;i++)
     {
        z=z+(gamma[i]-x[i])*V[i];
     }
     delete[] x;
     return {z};
}
