__kernel()
{
    int i;

    i=get_group_id(0);

    findinverse(Cloc[i].tensor);

    change_basis(xlambda6,stressG);
    change_basis(stressG6,stressG);
    change_basis(plasticStrain6,plasticStrain);
    getsymmetric(dispgrad);
    change_basis(strain6,strain);

//    sgnorm=0.
//    dgnorm=0.
//    do ii=1,3
//    do jj=1,3
//     sgnorm=sgnorm+xlambda(ii,jj)**2
//     dgnorm=dgnorm+strainaux(ii,jj)**2
//    enddo
//    enddo
//    sgnorm=sqrt(sgnorm)
//    dgnorm=sqrt(dgnorm)

    
}