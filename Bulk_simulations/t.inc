      parameter (maxnz=10001,maxNm=50001)
      real ran2
      common/tvec/xcp(maxNm),ycp(maxNm),zcp(maxNm),
     *xhp(maxNm),yhp(maxNm),zhp(maxNm),xcn(maxNm),ycn(maxNm),zcn(maxNm),
     *xhn(maxNm),yhn(maxNm),zhn(maxNm),rNsdist(maxNm),rNmdist(maxNm),
     *grhphp(maxnz),grhphn(maxnz),grhnhn(maxnz)
      common/ivec/iclustp(maxNm),iclustn(maxNm)
      common/var/ed2,ed,h,closew,q2,p2,utot,delta,
     *tx,ty,tz,cx,cy,cz,dpe,dpp,tex,tey,tez,edsq,ppf,psf,dpm,red2,
     *dmm2,Rsqcl,shuss,redsq,
     *ttx,tty,ttz,un,uo,red,bl,rbl,bl2,rbl2,
     *y5,y6,y7,y3,y4,Rcl,Rcl2,cldps,
     *rx,ry,rz,rrx,rry,rrz,qtx,qty,qtz,
     *dhspp,dhspp2,dhspn,dhspn2,dhsnn,dhsnn2,dpsp,dpsn,rnval,
     *rnval2,udn,udo,txs,tys,tzs,Rclust2,dpx,dpy,dpz,attclp,bjerrum,
     *rcut2,A4,u4,du4,rhspp,eratio,cdielc,dielc
      common/ihelt/Nm,Np,Nr,Ns,islu,ips,ipf,nhs,nhw,im,iref,ip,
     *idonk,ks,ms,kclp,kcln
