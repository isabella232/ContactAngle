/*-------------------------------------------------------------------------*\
You can redistribute this code and/or modify this code under the 
terms of the GNU General Public License (GPL) as published by the  
Free Software Foundation, either version 3 of the License, or (at 
your option) any later version. see <http://www.gnu.org/licenses/>.


The code has been developed by Ahmed AlRatrout as a part his PhD 
at Imperial College London, under the supervision of Dr. Branko Bijeljic 
and Prof. Martin Blunt. 

Please see our website for relavant literature:
AlRatrout et al, AWR, 2017 - https://www.sciencedirect.com/science/article/pii/S0309170817303342
AlRatrout et al, WRR, 2018 - https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2017WR022124
AlRatrout et al, PNAS, 2018 - http://www.pnas.org/

For further information please contact us by email:
Ahmed AlRarout:  a.alratrout14@imperial.ac.uk
Branko Bijeljic: b.bijeljic@imperial.ac.uk
Martin J Blunt:  m.blunt@imperial.ac.uk

Description
    Volume-preserving Gaussian smoothing

\*---------------------------------------------------------------------------*/


///. point-normal vectors
     pointField Cf=surf123.faceCentres();

     vectorField pNw(pointPoints.size(),vector(0,0,0));
     vectorField pNs(pointPoints.size(),vector(0,0,0));
     forAll(faces,faceI)
     {   const face & f=faces[faceI];
         forAll(f, pI)
         {   edge e=f.faceEdge(pI);
             vector Ce=0.5*(newPoints[e[0]] + newPoints[e[1]]); 

             vector pNE=0.5*((newPoints[e[0]]-Cf[faceI]) ^ (Ce-Cf[faceI])); 
             if( pMarks[e[0]]==fMarks[faceI] ) {    pNw[e[0]] += pNE;    pNs[e[0]] += pNE;    }
             else if( fMarks[faceI] != 2 )        pNw[e[0]] += pNE;
             else                                pNs[e[0]] += pNE;

             pNE=0.5*((newPoints[e[1]]-Cf[faceI]) ^ (Ce-Cf[faceI]));
             if( pMarks[e[1]]==fMarks[faceI] ) {    pNw[e[1]] -= pNE;    pNs[e[1]] -= pNE;    }
             else if( fMarks[faceI] != 2 )        pNw[e[1]] -= pNE;
             else                                pNs[e[1]] -= pNE;//2.*pNE-mag(pNE)*fNorms[faceI];
		 }    
     }
     scalarField pAreas(0.5*mag(pNs)+1e-18+0.5*mag(pNw));
     pNs/=(mag(pNs)+1e-18);
     pNw/=(mag(pNw)+1e-18);
		
	const vectorField previousPoints(newPoints);



	//===========================================================================
      Info<< iter   <<"  "; cout.flush();
      for (int iKern=0; iKern<kernelRadius; ++iKern)
      {
        vectorField displacements(newPoints.size(), vector(0,0,0));
        vectorField displacementsCL(newPoints.size(), vector(0,0,0));
        forAll(pointPoints, vertI)
        {

            const labelLoop& neiPoints = pointPoints[vertI];
            
			{
				vector avgPos(vector::zero);
				scalar sumWeights(1e-16); 
				if ( pMarks[vertI] == 4 )
				{
				  forAll(neiPoints, neiPointI)
				  {
					const label neiVertI = neiPoints[neiPointI];
					if( pMarks[neiVertI]!=2 )
					{
						//~ scalar weight=pWeights[neiVertI];//*max(NormalsOld[neiVertI]&NormalsOld[vertI],0.001);// /mag(newPoints[neiVertI]-newPoints[vertI]);
						scalar weight=(pAreas[neiVertI]+0.3*pAreas[vertI]);//*max(NormalsOld[neiVertI]&NormalsOld[vertI],0.001);// /mag(newPoints[neiVertI]-newPoints[vertI]);
						weight*=weight*pWeights[neiVertI];
						avgPos += weight * newPoints[neiVertI];
						sumWeights += weight;
					}
				  }

			  }
			  else
			  {
			    forAll(neiPoints, neiPointI)
				  {
					const label neiVertI = neiPoints[neiPointI];

					//~ scalar weight=pWeights[neiVertI];//*max(NormalsOld[neiVertI]&NormalsOld[vertI],0.001);// /mag(newPoints[neiVertI]-newPoints[vertI]);
					scalar weight=(pAreas[neiVertI]+0.3*pAreas[vertI]);//*max(NormalsOld[neiVertI]&NormalsOld[vertI],0.001);// /mag(newPoints[neiVertI]-newPoints[vertI]);
					weight*=weight*pWeights[neiVertI];
					avgPos += weight * newPoints[neiVertI];
					sumWeights += weight;

				  }
			  }
			  
				if (sumWeights>0.1)
				{
					avgPos /= sumWeights;//myEdges.size();
					displacements[vertI] = avgPos-newPoints[vertI]; //(0.8*relax)*((avgPos-newPoints[vertI])&pNw[vertI])*pNw[vertI] + (0.2*relax)*(avgPos-newPoints[vertI]);
				}
			}



            if( pMarks[vertI] == 4  )
            { ///. contact line smoothing
				vector avgPos(vector::zero);
				scalar sumWeights(1e-16); 
                forAll(neiPoints, pI)
                {
                    const label neiVertI = neiPoints[pI];
                    if( pMarks[neiVertI] == 4  )
                    {
                        scalar weight=1.;
                        avgPos += weight * newPoints[neiVertI];
                        sumWeights += weight;
                    }
                }

                vector DispCL =  (avgPos/sumWeights-newPoints[vertI]); // relax*((avgPos/sumWeights-newPoints[vertI])&pNs[vertI])*pNs[vertI]; 
                DispCL =  (DispCL & pNs[vertI])*pNs[vertI];
                displacementsCL[vertI] =  (DispCL - (DispCL & pNw[vertI])*pNw[vertI]);

            }

      }

		Info  <<"Gauss "; cout.flush();
		newPoints += (1.-extraRelaxFactorParIntrf)*relax*(displacements & pNw)*pNw + (extraRelaxFactorParIntrf*relax)*(displacements);
		newPoints += relaxCL*displacementsCL;

    }




	Info <<"     Volume-preservation "; cout.flush();

	vectorField dispAverage(newPoints-previousPoints);

  for (int iKern=0; iKern<kernelRadius+2; ++iKern)
  {
	const vectorField dispAverageTmp(dispAverage);

	//~ vectorField NormalsOld(pNw);
	forAll(pointPoints, vertI)
	{
		vector avgDisp(vector::zero);
		scalar sumWeights(1e-16);

		const labelLoop& neiPoints = pointPoints[vertI];
		label currentLabel = pMarks[vertI];
		bool curentLabelEq2 = currentLabel==2;

		forAll(neiPoints, neiPointI)
		{
			const label neiVertI = neiPoints[neiPointI];

			if( pMarks[neiVertI]==currentLabel || !( (pMarks[neiVertI]==2) ^ curentLabelEq2)  || pMarks[neiVertI] == 4  )
			{
				//~ scalar weight=1.;
				scalar weight=pAreas[neiVertI];
				if( pMarks[neiVertI]!=pMarks[vertI] ) weight*=0.05;
				avgDisp += weight * dispAverageTmp[neiVertI];
				sumWeights += weight;
			}
		}
		dispAverage[vertI] = avgDisp/sumWeights;//myEdges.size();
	  }
	}


	dispAverage -= 0.3*(dispAverage-(dispAverage& pNw)*pNw);
	dispAverage -= 0.3*(dispAverage-(dispAverage& pNs)*pNs);
	newPoints -= 1.09*dispAverage;     ///.  TODO change 0.5 0.5,  sum should be 1



	surf123.movePoints(newPoints);


	dispAverage=(newPoints-previousPoints);

	Info<< "  max:  "<<max(mag(dispAverage)); cout.flush();
	Info<< "  avgMag:  "<<average(mag(dispAverage)); cout.flush();
	Info<< "  avg:  "<<average(dispAverage); cout.flush();
	Info<<endl;

