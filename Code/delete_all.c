void delete_all()
{
    // These are data types defined in the "dirent" header
  /*  DIR* theFolder = opendir("./data");
    struct dirent *next_file;
    char filepath[256];
    */
    remove("./data/Pot.txt");
    remove("./data/Kin.txt");
    remove("./data/H.txt");
    remove("./data/vx.txt");
    remove("./data/vy.txt");
    remove("./data/vz.txt");
    remove("./data/rx.txt");
    remove("./data/ry.txt");
    remove("./data/rz.txt");
    remove("./data/r0x.txt");
    remove("./data/r1x.txt");
    remove("./data/v0x.txt");
    remove("./data/v1x.txt");
    remove("./data/gs.txt");
    remove("./data/h.txt");
    remove("./data/EigenvaluesLAPACK.txt");
    remove("./data/EigenvectorsLAPACK.txt");
    remove("./data/A.txt");
    remove("./data/Matrix.txt");
    remove("./data/alphaEnergy.txt");

}
