// Check for "points" in any of the result directories

bool hasMovingMesh = false;

if (timeDirs.size() > 1 && Pstream::master())
{
    // We already loaded a mesh (usually from constant).
    // See if any other "polyMesh/points" files exist too.

    Info<< "Search for moving mesh ... " << flush;
    for (const instant& inst : timeDirs)
    {
        const word& timeName = inst.name();

        hasMovingMesh =
        (
            timeName != mesh.pointsInstance()
         && IOobject
            (
                "points",
                timeName,
                polyMesh::meshSubDir,
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false  // no register
            ).typeHeaderOk<pointIOField>(true, false)
        );

        if (hasMovingMesh)
        {
            break;
        }
    }

    if (hasMovingMesh)
    {
        Info<< "found. Writing meshes for every timestep." << endl;
    }
    else
    {
        Info<< "none detected." << endl;
    }
}

reduce(hasMovingMesh, orOp<bool>());
