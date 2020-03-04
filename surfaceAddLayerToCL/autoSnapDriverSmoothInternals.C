#include "autoSnapDriver.H"
#include "motionSmoother.H"
#include "polyTopoChange.H"
#include "syncTools.H"
#include "fvMesh.H"
#include "Time.H"
#include "OFstream.H"
#include "mapPolyMesh.H"
#include "pointEdgePoint.H"
#include "PointEdgeWave.H"
#include "mergePoints.H"
#include "snapParameters.H"
#include "refinementSurfaces.H"
#include "unitConversion.H"

#include "primitivePatchInterpolation.H"

// Calculate displacement as average of patch points.
bool Foam::autoSnapDriver::smoothInternalPoints
(
    const snapParameters& snapParams,
    const label nInitErrors,
    motionSmoother& meshMover
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    pointVectorField& disp = meshMover.displacement();

    //~ const indirectPrimitivePatch& pp = meshMover.patch();
    //~ const polyMesh& mesh = meshMover.mesh();


    // Average cell centres
    // ~~~~~~~~~~~~~~


    const pointField& points = mesh.points();
    const vectorField& cellCentres = mesh.cellCentres();


{
    vectorField avgCellCentres(points.size(), vector::zero);
    labelList nPointsCells(points.size(), 0);

    forAll(cellCentres, cellI)
    {
        const labelList& 	cellPoints = mesh.cellPoints(cellI);

        forAll(cellPoints, pI)
        {
            label pointI = cellPoints[pI];
            avgCellCentres[pointI] += cellCentres[cellI];
            nPointsCells[pointI]++;
        }
    }

    syncTools::syncPointList
    (
        mesh,
        avgCellCentres,
        plusEqOp<point>(),  // combine op
        vector::zero        // null value
    );
    syncTools::syncPointList
    (
        mesh,
        nPointsCells,
        plusEqOp<label>(),  // combine op
        label(0)            // null value
    );

    forAll(points, i)
    {
        avgCellCentres[i] / nPointsCells[i];
        disp[i]=0.25*( avgCellCentres[i] / nPointsCells[i] - points[i] );
    }
}



//handle boundary faces;
/*{
    vectorField avgBFaceCentres(points.size(), vector::zero);
    labelList nBPointsFaces(points.size(), 0);    
    
    const fvBoundaryMesh& patches = mesh.boundary();
    forAll(patches, patchI)
    {
        if
        (
            ! patches[patchI].coupled()
        )
        {
            const polyPatch& pp=patches[patchI].patch();

            const vectorField::subField faceCentres = pp.faceCentres();

            forAll(pp, i)
            {
                const face& f = pp[i];
                const point& fc = faceCentres[i];

                forAll(f, fp)
                {
                    avgBFaceCentres[f[fp]] += fc;
                    nBPointsFaces[f[fp]]++;
                }
            }
        }
    }	 
    syncTools::syncPointList
    (
        mesh,
        avgBFaceCentres,
        plusEqOp<point>(),  // combine op
        vector::zero        // null value
    );
    syncTools::syncPointList
    (
        mesh,
        nBPointsFaces,
        plusEqOp<label>(),  // combine op
        label(0)            // null value
    );
    

     //~ avgBFaceCentres
     //~ nBPointsFaces
     
    forAll(patches, patchi)
    {
        const polyPatch& pp =patches[patchi].patch();
        vectorField pPointNormals(patches[patchi].patch().pointNormals());
        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),           
            pPointNormals,
            plusEqOp<point>(),  // combine op
            vector::zero        // null value
        );
        pPointNormals/=mag(pPointNormals)+1e-12;
		//~ primitivePatchInterpolation pinterpolator(pp);
        //~ vectorField facesCentres=pp.faceCentres();
        //~ vectorField faceCentresAvg=pinterpolator.faceToPointInterpolate(facesCentres);
        forAll(pPointNormals, pi)
        {
            label pointI=pp.meshPoints()[pi];
            if(nBPointsFaces[pointI]>1)
            {
                disp[pointI] = 0.5*(pPointNormals[pi]^(disp[pointI]^pPointNormals[pi]));
                                        //~ +0.1*(faceCentresAvg[pi]-points[pointI]);               
            }
            else
            {
                disp[pointI] = 0.5*(pPointNormals[pi]^((avgBFaceCentres[pi]-points[pointI])^pPointNormals[pi]));
                            
            }
        }
   
    }	 
    



}*/




    //~ forAll(points,i)
    //~ {
        //~ meshMover.displacement()[i]=0.5*(avgCellCentres[i]-points[i]);
    //~ }
    meshMover.displacement().boundaryField()==vector(0.0,0.0,0.0);

    meshMover.displacement().boundaryField().evaluate();
    meshMover.displacement().correctBoundaryConditions();
    

    
    List<labelPair> emptyBaffles;
    bool meshOK = scaleMesh
    (
        snapParams,
        nInitErrors,
        emptyBaffles,
        meshMover
    );
 
    //~ mesh.movePoints(avgCellCentres);

    return meshOK;
} 










/*


void Foam::setDisplacement(pointField& displacement)
{
    // See comment in .H file about shared points.
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            const labelList& meshPoints = pp.meshPoints();

            forAll(meshPoints, i)
            {
                displacement_[meshPoints[i]] = vector::zero;
            }
        }
    }

    const labelList& ppMeshPoints = pp_.meshPoints();

    // Set internal point data from displacement on combined patch points.
    forAll(ppMeshPoints, patchPointI)
    {
        displacement_[ppMeshPoints[patchPointI]] = displacement[patchPointI];
    }

    // Adapt the fixedValue bc's (i.e. copy internal point data to
    // boundaryField for all affected patches)
    forAll(adaptPatchIDs_, i)
    {
        label patchI = adaptPatchIDs_[i];

        displacement_.boundaryField()[patchI] ==
            displacement_.boundaryField()[patchI].patchInternalField();
    }

    // Make consistent with non-adapted bc's by evaluating those now and
    // resetting the displacement from the values.
    // Note that we're just doing a correctBoundaryConditions with
    // fixedValue bc's first.
    labelHashSet adaptPatchSet(adaptPatchIDs_);

    const lduSchedule& patchSchedule = mesh_.globalData().patchSchedule();

    forAll(patchSchedule, patchEvalI)
    {
        label patchI = patchSchedule[patchEvalI].patch;

        if (!adaptPatchSet.found(patchI))
        {
            if (patchSchedule[patchEvalI].init)
            {
                displacement_.boundaryField()[patchI]
                    .initEvaluate(Pstream::scheduled);
            }
            else
            {
                displacement_.boundaryField()[patchI]
                    .evaluate(Pstream::scheduled);
            }
        }
    }

    // Multi-patch constraints
    applyCornerConstraints(displacement_);

    // Correct for problems introduced by corner constraints
    syncTools::syncPointList
    (
        mesh_,
        displacement_,
        maxMagEqOp(),   // combine op
        vector::zero    // null value
    );

    // Adapt the fixedValue bc's (i.e. copy internal point data to
    // boundaryField for all affected patches) to take the changes caused
    // by multi-corner constraints into account.
    forAll(adaptPatchIDs_, i)
    {
        label patchI = adaptPatchIDs_[i];

        displacement_.boundaryField()[patchI] ==
            displacement_.boundaryField()[patchI].patchInternalField();
    }

    if (debug)
    {
        OFstream str(mesh_.db().path()/"changedPoints.obj");
        label nVerts = 0;
        forAll(ppMeshPoints, patchPointI)
        {
            const vector& newDisp = displacement_[ppMeshPoints[patchPointI]];

            if (mag(newDisp-displacement[patchPointI]) > SMALL)
            {
                const point& pt = mesh_.points()[ppMeshPoints[patchPointI]];

                meshTools::writeOBJ(str, pt);
                nVerts++;
                //Pout<< "Point:" << pt
                //    << " oldDisp:" << displacement[patchPointI]
                //    << " newDisp:" << newDisp << endl;
            }
        }
        Pout<< "Written " << nVerts << " points that are changed to file "
            << str.name() << endl;
    }

    // Now reset input displacement
    forAll(ppMeshPoints, patchPointI)
    {
        displacement[patchPointI] = displacement_[ppMeshPoints[patchPointI]];
    }
}



*/









//~ 
//~ 
//~ 
//~ void Foam::autoSnapDriver::doSmoothInternalPoints
//~ (
    //~ const dictionary& snapDict,
    //~ const dictionary& motionDict,
    //~ const scalar featureCos,
    //~ const snapParameters& snapParams
//~ )
//~ {
    //~ fvMesh& mesh = meshRefiner_.mesh();
//~ 
    //~ Info<< nl
        //~ << "Morphing phase" << nl
        //~ << "--------------" << nl
        //~ << endl;
//~ 
    //~ // Get the labels of added patches.
    //~ labelList adaptPatchIDs(meshRefiner_.meshedPatches());
//~ 
    //~ // Create baffles (pairs of faces that share the same points)
    //~ // Baffles stored as owner and neighbour face that have been created.
    //~ List<labelPair> baffles;
    //~ meshRefiner_.createZoneBaffles(globalToPatch_, baffles);
//~ 
//~ 
    //~ bool doFeatures = false;
    //~ label nFeatIter = 1;
    //~ if (snapParams.nFeatureSnap() > 0)
    //~ {
        //~ doFeatures = true;
        //~ nFeatIter = snapParams.nFeatureSnap();
//~ 
        //~ Info<< "Snapping to features in " << nFeatIter
            //~ << " iterations ..." << endl;
    //~ }
//~ 
//~ 
    //~ bool meshOk = false;
//~ 
    //~ {
        //~ autoPtr<indirectPrimitivePatch> ppPtr
        //~ (
            //~ meshRefinement::makePatch
            //~ (
                //~ mesh,
                //~ adaptPatchIDs
            //~ )
        //~ );
//~ 
        //~ // Distance to attract to nearest feature on surface
        //~ const scalarField snapDist(calcSnapDistance(snapParams, ppPtr()));
//~ 
//~ 
        //~ // Construct iterative mesh mover.
        //~ Info<< "Constructing mesh displacer ..." << endl;
        //~ Info<< "Using mesh parameters " << motionDict << nl << endl;
//~ 
        //~ const pointMesh& pMesh = pointMesh::New(mesh);
//~ 
        //~ motionSmoother meshMover
        //~ (
            //~ mesh,
            //~ ppPtr(),
            //~ adaptPatchIDs,
            //~ meshRefinement::makeDisplacementField(pMesh, adaptPatchIDs),
            //~ motionDict
        //~ );
//~ 
//~ 
//~ 
//~ 
//~ 
//~ 
     //~ Info<< "\nV''''''''''''''''''''''''''''''V''''''''''''''''''''''''''''''V\n"<<0<<endl;
        //~ smoothInternalPoints            
            //~ (
                //~ snapParams,
                //~ 1000000000,
                //~ meshMover
            //~ );
    //~ Info<< "A_____________________________A________________________________A\n";
            //~ meshMover.correct();
//~ 
//~ 
     //~ Info<< "\nV''''''''''''''''''''''''''''''V''''''''''''''''''''''''''''''V\n"<<0<<endl;
        //~ smoothInternalPoints            
            //~ (
                //~ snapParams,
                //~ 1000000000,
                //~ meshMover
            //~ );
    //~ Info<< "A_____________________________A________________________________A\n";
            //~ meshMover.correct();
//~ 
    //~ 
    //~ 
    //~ 
//~ 
        //~ // Check initial mesh
        //~ Info<< "Checking initial mesh ..." << endl;
        //~ labelHashSet wrongFaces(mesh.nFaces()/100);
        //~ motionSmoother::checkMesh(false, mesh, motionDict, wrongFaces);
        //~ const label nInitErrors = returnReduce
        //~ (
            //~ wrongFaces.size(),
            //~ sumOp<label>()
        //~ );
//~ 
        //~ Info<< "Detected " << nInitErrors << " illegal faces"
            //~ << " (concave, zero area or negative cell pyramid volume)"
            //~ << endl;
//~ 
//~ 
        //~ Info<< "Checked initial mesh in = "
            //~ << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;
//~ 
 //~ 
 //~ 
     //~ Info<< "\nV''''''''''''''''''''''''''''''V''''''''''''''''''''''''''''''V\n"<<nInitErrors<<endl;
        //~ smoothInternalPoints            
            //~ (
                //~ snapParams,
                //~ nInitErrors,
                //~ meshMover
            //~ );
    //~ Info<< "A_____________________________A________________________________A\n";
             //~ meshMover.correct();
//~ 
 //~ 
        //~ // Pre-smooth patch vertices (so before determining nearest)
        //~ preSmoothPatch(snapParams, nInitErrors, baffles, meshMover);
//~ 
    //~ Info<< "\nV''''''''''''''''''''''''''''''V''''''''''''''''''''''''''''''V\n"<<nInitErrors<<endl;
        //~ smoothInternalPoints            
            //~ (
                //~ snapParams,
                //~ nInitErrors,
                //~ meshMover
            //~ );
    //~ Info<< "A_____________________________A________________________________A\n";
            //~ meshMover.correct();
//~ 
        //~ for (label iter = 0; iter < nFeatIter; iter++)
        //~ {
        //~ }
    //~ }
//~ 
//~ }













