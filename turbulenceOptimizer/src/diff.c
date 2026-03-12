
scalar diff
(
  volScalarField &gamma,
  const scalarField &V,
  double del,
  double eta,
  int n,
  const Switch projectionNumericalSafeguards
)
{
     int i;
     scalar z=0;
     double *x =new double[n];
     const scalar etaEps(1e-6);
     const scalar etaProjection =
       projectionNumericalSafeguards
       ? Foam::min(Foam::max(eta, etaEps), scalar(1.0)-etaEps)
       : eta;
     const scalar etaDenominator =
       projectionNumericalSafeguards ? Foam::max(etaProjection, etaEps) : etaProjection;
     const scalar oneMinusEtaDenominator =
       projectionNumericalSafeguards
       ? Foam::max(scalar(1.0)-etaProjection, etaEps)
       : (scalar(1.0)-etaProjection);
     
     for(i=0;i<n;i++)
     {
        const scalar gammaProjection =
          projectionNumericalSafeguards
          ? Foam::min(Foam::max(gamma[i], scalar(0.0)), scalar(1.0))
          : gamma[i];
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
