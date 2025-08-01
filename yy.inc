      parameter (maxk=30001,maxnz=10001,
     *maxNm=30001,kkmax=300001)
      real ran2
      common/vec/xm(maxNm),ym(maxNm),zm(maxNm),
     *wp(maxNm),ttwp(maxNm),
     *xs(maxnm),ys(maxnm),zs(maxnm),wps(maxnm),ttwps(maxNm),
     *fvec(kkmax),
     *vzrk(kkmax),vyrk(kkmax),vxrk(kkmax),
     *vzcsum(kkmax),vzssum(kkmax),vcsum(kkmax),
     *vssum(kkmax),
     *tvcsum(kkmax),tvssum(kkmax),
     *tfvec(kkmax),
     *tvxrk(kkmax),tvzrk(kkmax),
     *grcc(maxnz),grca(maxnz),graa(maxnz),
     *xcn(maxnm),ycn(maxnm),zcn(maxnm),
     *xcp(maxnm),ycp(maxnm),zcp(maxnm)

      common/ivec/iclustp(maxNm),iclustn(maxNm)

      common/var/ed2,ed,h,red,umm,uw,rNm,rNp,rNr,rNs,
     *uwm,awm,rwm,ttuw,umn,umo,duw,tx,ty,tz,rnn,
     *tex,tey,tez,cx,cy,cz,dpee,dpe,dpp,red2,chemp,twopi,
     *volume,Wtot,rktrial,ene,eno,eps,Q,epsmod,sigma,
     *alpha,a1,a2,a3,a4,a5,p,vk,vr,ra,fup,p1,p2,p3,
     *ur,uk,deltauk,deltaur,pi,fourpi,efact,rkfact,bjerrum,
     *edsq,rvolume,dmm2,dss2,dms2,rnmval,rnsval,rnmval2,
     *dps,halfh,rh,rkzfact,urself,bbur,rtz,rttz,durself,
     *zrw,zlw,uab,deltauab,quarth,umi,dumi,ubbmi,dubbmi,
     *ukn,urn,rhs,rhs2,rrhs,pzrw,pzlw,fzk,fzr,
     *fzlhsr,fzrhsr,fzlhsk,fzrhsk,urljn,urljo,urlj,
     *frljl,frljr,dUrlj,dUr,dUk,trkzfact,tuk,tur,trvolume,
     *udz,udh,ttur,ttuk,th,trh,turlj,tturlj,fracslump,
     *Rclust2,dpx,dpy,dpz,txs,tys,tzs,rnval,rnval2,
     *udn,udo,cldps,rhspp,eratio,cdielc,dielc,dhspp,rcut
    
      common/iheltal/Nm,Np,Nr,idum,ks,ip,im,ips,ipf,
     *idonk,ktrial,l,nn,n,ns,iz,j,kqq,knt0,kntmax,
     *kmax,ksqmax,kntotmax,kzmax,nrefl,ms,kclp,kcln  



