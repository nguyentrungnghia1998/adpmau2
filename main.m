ug1=0;
ug2=0;
ug0=0;
yg0=1;
yg1=1;
yg2=1;
Rd=1;
Qd=1;
gamma=0.2;
siold=[3;4];
Pd_=ones(4);
ukdata=zeros(100,1);
ykdata=zeros(100,1);
for k=1:100
    upre1=ug1;
    ypre=[yg0;yg1];
    ek=nhieu(k);
    ug0=-1/(Rd/gamma+Pd_(1,1))*(Pd_(1,2)*upre1+Pd_(1,3:4)*ypre)+ek;
    yg2=yg1;
    yg1=yg0;
    [yg0,sinew]=mohinh(siold,ug0);
    ug2=ug1;
    ug1=ug0;
    siold=sinew;
    ukdata(k)=ug0;
    ykdata(k)=yg0;
end
Odata=zeros(98,16);
for k=3:100
    upre=ukdata(k-1:-1:k-2);
    ypre=ykdata(k-1:-1:k-2);
    zpre=[upre;ypre];
    Odata(k-2,:)=kron(zpre',zpre');
end
ujdata=ukdata;
ujdatanext=zeros(98,1);
for j=0:100
    Vdata=zeros(98,1);
    for k=3:100
        unext=ukdata(k:-1:k-1);
        ynext=ykdata(k:-1:k-1);
        znext=[unext;ynext];
        Vdata(k-2)=ykdata(k)'*Qd*ykdata(k)+ujdata(k)'*Rd*ujdata(k)+gamma*znext'*Pd_*znext;
    end
    Pnewstk=pinv(Odata)*Vdata;
    Pnew=reshape(Pnewstk,[4 4]);
    if norm(Pnew-Pd_,1)<0.0005
        break;
    end
    for k=3:100
        upre1=ukdata(k-1);
        ynext=ykdata(k:-1:k-1);
        ujdata(k-2)=-1/(Rd/gamma+Pnew(1,1))*(Pnew(1,2)*upre1+Pnew(1,3:4)*ynext);
    end
    Pd_=Pnew
end