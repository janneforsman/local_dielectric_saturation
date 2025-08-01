      implicit double precision (a-h,o-z)      
      include 'yy.inc'              
      dimension rhozm(maxnz),rhozs(maxnz),rNsdist(maxNm),rNmdist(maxNm) 
      pi = acos(-1.d0)
      twopi = 2.d0*pi
      fourpi = 4.d0*pi
      dpe = twopi
      dpee = twopi
      rtwelve = 1.d0/12.d0
      rthree = 1.d0/3.d0
      volfact = fourpi/3.d0
      rnine = 1.d0/9.d0
      rthree = 1.d0/3.d0
      elch = 1.602D-19
      avno = 6.02214D23
      bk = 1.38066D-23
      faraday = 8.85418782D-12
      temp = 298.d0
c      temp = 561.d0

      a1 = 0.254829592d0
      a2 = -0.284496736d0
      a3 = 1.421413741d0
      a4 = -1.453152027d0
      a5 = 1.061405429d0      

      do i = 1,maxnz
      rhozs(i) = 0.d0
      rhozm(i) = 0.d0   
      grcc(i) = 0.d0
      grca(i) = 0.d0
      graa(i) = 0.d0
      enddo
      maxNp = 1
      maxNn = 1
      minNn = maxNm
      minNp = maxNm      
      do i = 1,maxNm
      rNsdist(i) = 0.   
      rNmdist(i) = 0.   
      enddo      
      Nnd = 11
      Npd = 12      
      kkk = 85
      lll = 87   
      mmm = 89   
      islump = 96
      irdf = 92
      inzp = 54
      inzn = 56
      inf = 29
      open(islump,file='slump',form='formatted')
      rewind islump
      read(islump,*) idum
      write(*,*) 'idum,ran2(idum) = ',idum,ran2(idum)
      write(*,*) 'idum,ran2(idum) = ',idum,ran2(idum)
      open (inf,file='Nfil',form='formatted')        
      open(mmm, file='input',form='formatted')
      open(lll, file='coordm',form='formatted')
      open(kkk, file='coords',form='formatted')
      open(irdf, file='rdf',form='formatted')
      open(inzp, file='nzm',form='formatted')
      open(inzn, file='nzs',form='formatted')
      open (Nnd,file='Nsdist',form='formatted')
      open (Npd,file='Nmdist',form='formatted')      
      rewind mmm
      read(mmm,*) 
      read(mmm,*) Nsub
      read(mmm,*) 
      read(mmm,*) Nmtrams,nmvaltrams
      read(mmm,*) 
      read(mmm,*) Nstrams
      write(*,*) Nstrams
      nsval = 1
c     (s is implcitly assumed to be negative, and monovalent)
      read(mmm,*)
c     pdw = closest distance of approach for cations (0.5*dhs for anions)      
      read(mmm,*) ed2,h
      write(*,*) ed2,h
c     h = box length along z, ed2 = box length along x,y
      read(mmm,*) 
      read(mmm,*) dpp,dps
      read(mmm,*) 
      read(mmm,*) deltaz
      read(mmm,*) 
      read(mmm,*) knsamp,nspot
      read(mmm,*) 
      read(mmm,*) kread
      read(mmm,*) 
      read(mmm,*) dmm,dss,dms
      dss = dmm
      dms = dmm
      write(*,*) 'dss = dms = dmm = ',dmm
      read(mmm,*) 
      read(mmm,*) dielc,temp,cdielc 
      write(*,*) 'dielc,temp = ',dielc,temp
      write(*,*) 'dielc,cdielc = ',dielc,cdielc
      diffdielc = dielc-cdielc      
      read(mmm,*) 
      read(mmm,*) fracmval,fracgcmove,fracslump
      read(mmm,*) 
      read(mmm,*) alphapar,kmax,kzmax
      read(mmm,*) 
      read(mmm,*) fracclm,Rclmax,cldps
      write(*,*) 'OBS.! Varying cluster radius, Rclmax = ',Rclmax          
      write(*,*) 'cldps = ',cldps
      write(*,*) 'fracclm = ',fracclm                  

      rewind inf
      read(inf,*) Ns
      read(inf,*) Nm,nmval
      read(inf,*) chemp
      read(inf,*) donnan      
      
c      ksqmax = 3*kmax**2
c     kntotmax=(2*kmax+1)**3
      ksqmax = 2*kmax**2      
c     kntotmax=(kmax+kzmax+1)**3
      kntotmax = (2*kmax+1)*(2*kzmax+1)*(kmax+1)
      nchelec = Nm*nmval-Ns*nsval
c      if (nchelec.ne.0) goto 9999
      rnsval = dfloat(Nsval)
      rnmval = dfloat(Nmval)      
      rnmval2 = rnmval**2
      rnval = rnmval
      rnval2 = rnmval2
      sigmamet = 1.D-10
      rtemp = dielc*4.d0*pi*faraday*sigmamet*bk*temp/elch**2
      rbeta = 1.d0/rtemp
      write(*,*) 'rtemp,rbeta = ',rtemp,rbeta
      bjerrum = rbeta
      dmm2 = dmm*dmm
      dss2 = dss*dss
      dms2 = dms*dms

      dhspp = dmm
      rcut = 2.d0*dhspp      
      rcut2 = rcut**2
      rhspp = dhspp/2.d0      
      
      volume  = ed2*ed2*h
      rvolume = 1.d0/volume
      ed = 0.5d0*ed2
      halfh = 0.5d0*h
      rh = 1.d0/h
      red = 1.d0/ed
      red2 = 1.d0/ed2
      edsq = ed*ed
      rarea = 1.d0/ed2**2
      subvolume = ed2*ed2*h
      psubvolume = ed2*ed2*h
      bdm = real(Nm)/psubvolume
      bds = real(Ns)/subvolume      
      write(*,*) '********************************'
      write(*,*) 'Simulation of bulk salt (EWALD)'
      write(*,*) '********************************'
      write(*,*) 'no. of configurations/1000: ',Nsub/100
      write(*,*) 'Nm,Ns = ',Nm,Ns
      write(*,*) 'Nmval,Nsval = ',Nmval,Nsval
      write(*,*) 'dmm,dss = ',dmm,dss
      write(*,*) 'dms = ',dms
      write(*,*) 'dielc = ',dielc
      write(*,*) 'rcut = dhspp = ',rcut,rcut2      
      write(*,*) 'epsr = cdielc+diffdielc*(r-dhspp)/dhspp, r.lt.rcut'  
      eratio = diffdielc/dhspp
      write(*,*) 'eratio,diffdielc = ',eratio,diffdielc        
      write(*,*) 'fracmval = ',fracmval
      write(*,*) 'fracgcmove = ',fracgcmove
      write(*,*) 'fracslump = ',fracslump    
      write(*,*) 'lateral length of box along x,y (ed2): ',ed2
      write(*,*) 'lateral length of box along z (h): ',h      
      write(*,*) 'dpp,dps: ',dpp,dps
      write(*,*) 'knsamp  = ',knsamp
      write(*,*) 'dr (g(r) distr.) = ',deltaz      
      alpha = alphapar/ed2
      write(*,*) 'alphapar = ',alphapar
      write(*,*) 'alpha = ',alpha
      write(*,*) 'kmax,kzmax = ',kmax,kzmax
      write(*,*) 'ksqmax = ',ksqmax
      write(*,*) 'kntotmax = ',kntotmax      
      if (kread.eq.1) then
      rewind lll
      do  i = 1,Nm
      read(lll,*) xm(i),ym(i),zm(i)
      enddo
      rewind kkk
      do i = 1,Ns
      read(kkk,*) xs(i),ys(i),zs(i)
      enddo      
      elseif (kread.eq.2) then           
      do i = 1,Nm
      read(kkk,*) xm(i),ym(i),zm(i)
      read(kkk,*) xs(i),ys(i),zs(i)
      enddo
      else
      call cinit
      endif
      rNs = dfloat(Ns)
      rNm = dfloat(Nm)      
      efact = 0.25d0/alpha**2
      rkfact = 2.d0*pi/ed2
      rkzfact = 2.d0*pi/h
      call ksetup
      call kvecs
      call check
c     twice as many ions:
      sfactor = (rNm*rnmval**2+rNs*rnsval**2)
      uself = -sfactor*rbeta*alpha/dsqrt(pi)
      write(*,*) 'uk,ur = ',uk,ur
      write(*,*) 'umi = ',umi
      write(*,*) 'uself = ',uself
      write(*,*) 'uk+ur= ',uk+ur
c      write(*,*) 'urlj = ',urlj
      write(*,*) 'utot = ',uk+ur+uself
      write(*,*) 
      write(*,*) 'salt chemical potential: ',chemp
      write(*,*) 'Donnan potential*e:',Donnan
      chs = chemp+Donnan
      chm = chemp-Donnan*rnmval
      write(*,*) 'chm,chs = ',chm,chs
      write(*,*) 'fracgcmove = ',fracgcmove
      write(*,*) 'rnmval,rnsval = ',rnmval,rnsval      

c     when adding an m ion:
      duselfm = -rnmval**2*rbeta*alpha/dsqrt(pi)
c     when adding an s ion:      
      duselfs = -rnsval**2*rbeta*alpha/dsqrt(pi)            
      
      Nmac = 10
      rdeltaz = 1.d0/deltaz
      dr = deltaz
      rdr = rdeltaz
      rewind 98
      r = -0.5d0*deltaz
      ndr = nint(min(ed,h)/deltaz)
      do i = 1,ndr
      r = r+deltaz
      if (r.gt.dmm) write(98,*) r,terfc(alpha*r)/r,erfc(alpha*r)/r
      enddo

      rewind 78
      r = dhspp-dr
      nnn = nint(10./dr)+1
      do i = 1,nnn
      r = r+dr
      rmm2 = r**2
      if (rmm2.lt.rcut2) then
      r = dsqrt(rmm2)
      tdielc = cdielc+eratio*(r-dhspp)      
      ppot = dielc/tdielc-1.d0
c      udo = -rnval*(dielc/tdielc-1.d0)+udo      
      else
      ppot = 0.d0
      endif               
      pot4 = bjerrum*ppot/r      
      write(78,*) r,pot4,pot4+bjerrum/r
      enddo
c      stop
      
      asampnz = 0.d0
      accs = 0.
      rejs = 0.
      atts = 0.
      attm = 0.
      accm = 0.
      rejm = 0.
      rnms = 0.
      rnmm = 0.
      rnss = 0.
      tvfr = 4.d0*pi/3.d0
      vsph = (4.d0*pi*ed**3)/3.d0

      attcrem = 0.
      acccrem = 0.
      rejcrem = 0.
      attcres = 0.
      acccres = 0.
      rejcres = 0.
      attdesm = 0.
      accdesm = 0.
      rejdesm = 0.
      attdess = 0.
      accdess = 0.
      rejdess = 0.
      avNs = 0.
      avNm = 0.

      surfpot = 0.d0
      cpot = 0.d0
      epsilon = dielc
      epszero = faraday
      dhs = dmm
      cnorm = elch/(epsilon*epszero*1.D-10)      
      pnorm = cnorm/fourpi
      rhs = 0.5d0*dhs
      rrhs = 1.d0/rhs
      rhs2 = rhs**2
      write(*,*) 'cnorm,pnorm = ',cnorm,pnorm
      write(*,*) 'rhs,nspot = ',rhs,nspot
      write(*,*) 'bk*temp/elch = ', bk*temp/elch      
c      call flush
      sumcpot = 0.d0
      sumspot = 0.d0
      savesampsdens = 0.d0
      saveavdens = 0.d0
      avsdens = 0.d0
      sampsdens = 0.d0
      savesdens = 0.d0
      samper = 0.d0
      sduk = 0.d0
      sdurlj = 0.d0
      sdur = 0.d0
      sedVim = 0.d0
c      write(*,*) 'deltah?'
c      read(*,*) deltah
      deltah = 0.001d0
      write(*,*) 'deltah = ',deltah
      th = h-deltah
c      th = h-2.d0*deltah
      write(*,*) 'th = ',th
      trkzfact =  2.d0*pi/th
      trh = 1.d0/th
      tvolume = th*ed2**2
      trvolume = 1.d0/tvolume
      halfdh = 0.5d0*deltah
      dV = ed2*ed2*deltaz
      ihalfh = nint(halfh/deltaz)
      sdens = 0.5d0*red2**2*(rNs-rnmval*rNm)      
      attclp = 0.
      attcln = 0.
      accclp = 0.
      acccln = 0.
      rejclp = 0.
      rejcln = 0.      
      do 5 kmac = 1,Nmac
      do 276 ksub = 1,Nsub        
      sampsdens = 1.+sampsdens
      avsdens = sdens+avsdens
      idonk = 0
      slumpt = ran2(idum)
      
      if (slumpt.lt.fracgcmove) then
c     GC step:
      if (ran2(idum).lt.0.5d0) then
c     attempted GC step for an m ion
         
      if (ran2(idum).lt.0.5d0) then
c     attempted addition of an m ion
      attcrem = attcrem+1.
      call crem
      
      if (idonk.eq.1) goto 7687
      ukn = deltauk
      urn = deltaur
      delta = ukn+urn+duselfm
      ugc = delta
      exponent = (chm-ugc)+dlog(psubvolume/(rNm+1.d0))      
      if (exponent.lt.-120.d0) goto 7687
      prob = dexp(exponent)      
      if (prob.lt.ran2(idum)) goto 7687
c     accepted creation
      acccrem = acccrem+1.
      Nm = Nm+1
      xm(Nm) = tx
      ym(Nm) = ty
      zm(Nm) = tz
      uk = uk+ukn
      ur = ur+urn
      uself = uself+duselfm
c      urlj = urlj+urljn
      rNm = dfloat(Nm)
      do kntot = 1,kntmax
      vcsum(kntot) = tvcsum(kntot)
      vssum(kntot) = tvssum(kntot)
      enddo           
      goto 6
 7687 rejcrem = rejcrem+1.
      goto 6

      else
c     attempted destruction of an m ion
      attdesm = attdesm+1
      if (Nm.eq.1) goto 7689
      ks = int(1.d0+rNm*ran2(idum))
      if (ks.gt.Nm) ks = Nm
      if (ks.lt.1) ks = 1         
      call desm
      ukn = -deltauk      
      urn = deltaur
      delta = ukn+urn+duselfm
      ugc = delta
c     ugc is interaction energy of removing particle with surroundings
      exponent = (-chm+ugc)+dlog(rNm/psubvolume)
      if (exponent.lt.-120.d0) goto 7689
      prob = dexp(exponent)            

      if (prob.lt.ran2(idum)) goto 7689
c     accepted destruction
      accdesm = accdesm+1.
      uk = uk-ukn
      ur = ur-urn
      uself = uself-duselfm
c      urlj = urlj-urljn
      Nm = Nm-1
      do iii = ks,Nm
      xm(iii) = xm(iii+1)
      ym(iii) = ym(iii+1)
      zm(iii) = zm(iii+1)
      enddo
      rNm = dfloat(Nm)
      do kntot = 1,kntmax
      vcsum(kntot) = tvcsum(kntot)
      vssum(kntot) = tvssum(kntot)
      enddo           
      goto 6
 7689 rejdesm = rejdesm+1.
      goto 6
      endif

      else
c     attempted GC step for an s ion
         
      if (ran2(idum).lt.0.5d0) then
c     attempted addition of an n ion
      attcres = attcres+1.      
      call cres
      if (idonk.eq.1) goto 7627

      ukn = deltauk
      urn = deltaur
      delta = ukn+urn+duselfs
      ugc = delta
      exponent = (chs-ugc)+dlog(subvolume/(rNs+1.d0))
      if (exponent.lt.-120.d0) goto 7627
      prob = dexp(exponent)      
      if (prob.lt.ran2(idum)) goto 7627
c     accepted creation
      acccres = acccres+1.
      Ns = Ns+1
      xs(Ns) = tx
      ys(Ns) = ty
      zs(Ns) = tz
      uk = uk+ukn
      ur = ur+urn      
      uself = uself+duselfs
c      urlj = urlj+urljn
      rNs = dfloat(Ns)
      do kntot = 1,kntmax
      vcsum(kntot) = tvcsum(kntot)
      vssum(kntot) = tvssum(kntot)
      enddo      
      goto 6
 7627 rejcres = rejcres+1.
      goto 6

      else
c     attempted destruction of an n ion
      attdess = attdess+1
      if (Ns.eq.1) goto 7629
      ks = int(1.d0+rNs*ran2(idum))
      if (ks.gt.Ns) ks = Ns
      if (ks.lt.1) ks = 1         
      call dess
      ukn = -deltauk
      urn = deltaur
      delta = ukn+urn+duselfs
      ugc = delta
c     ugc is interaction energy of removing particle with surroundings
      exponent = (-chs+ugc)+dlog(rNs/subvolume)
      if (exponent.lt.-120.d0) goto 7629
      prob = dexp(exponent)            

      if (prob.lt.ran2(idum)) goto 7629
c     accepted destruction
      accdess = accdess+1.
      uk = uk-ukn
      ur = ur-urn
      uself = uself-duselfs
c      urlj = urlj-urljn
      Ns = Ns-1
      do iii = ks,Ns
      xs(iii) = xs(iii+1)
      ys(iii) = ys(iii+1)
      zs(iii) = zs(iii+1)
      enddo
      rNs = dfloat(Ns)
      sdens = 0.5d0*red2**2*(rNs-rnmval*rNm)
      do kntot = 1,kntmax
      vcsum(kntot) = tvcsum(kntot)
      vssum(kntot) = tvssum(kntot)
      enddo
      
      goto 6
 7629 rejdess = rejdess+1.
      goto 6
      endif
      endif
      
      else
c     normal displacements:
      
      if (ran2(idum).lt.fracmval) then
c     attempted move of m ion:

      if (ran2(idum).lt.fracclm) then
c     attempted cluster displacement of cation
      ks = int(1.d0+rNm*ran2(idum))
      if (ks.gt.Nm) ks = Nm
      if (ks.lt.1) ks = 1                           
      Rclust = dmm+ran2(idum)*(Rclmax-dmm)
      Rclust2 = Rclust*Rclust
      attclp = attclp+1.
c      write(*,*) 'attclp = ',attclp
      call clp
      if (idonk.eq.1) goto 927
      dut = deltaur+deltauk
      if (dut.gt.100.d0) goto 927      
      if (dut.lt.0.d0) goto 9165
      if (dexp(-dut).lt.ran2(idum)) goto 927
 9165 accclp = accclp+1.      
      ur = ur+deltaur
      xm(ks) = txs
      ym(ks) = tys
      zm(ks) = tzs
      do kk = 1,kcln
      km = iclustn(kk)
      xs(km) = xcn(km)
      ys(km) = ycn(km)
      zs(km) = zcn(km)
      enddo
      do kk = 1,kclp
      km = iclustp(kk)
      xm(km) = xcp(km)
      ym(km) = ycp(km)
      zm(km) = zcp(km)
      enddo
      uk = uk+deltauk
      do kntot = 1,kntmax
      vcsum(kntot) = tvcsum(kntot)
      vssum(kntot) = tvssum(kntot)
      enddo      
      goto 6
 927  rejclp = rejclp+1.
      goto 6
      
      else
c     attempted standard displacement of cation:      
      attm = attm+1.
      call mtrans
      if (idonk.eq.1) goto 103
c      deltauk = 0.5d0*deltauk
c      deltaur = 0.5d0*deltaur
c      deltaurlj = urljn-urljo      
      deltau = deltauk+deltaur
      if (dexp(-deltau).lt.ran2(idum)) goto 103
      accm = accm+1.
      uk = uk+deltauk
      ur = ur+deltaur
c      urlj = urlj+deltaurlj
      xm(ip) = tx
      ym(ip) = ty
      zm(ip) = tz
      do kntot = 1,kntmax
      vcsum(kntot) = tvcsum(kntot)
      vssum(kntot) = tvssum(kntot)
      enddo
      
      goto 105
 103  rejm = rejm+1.
 105  goto 6
      endif
      
      else
c     attempted move of s ion:

      if (ran2(idum).lt.fracclm) then
c     attempted cluster displacement of cation
      ks = int(1.d0+rNm*ran2(idum))
      if (ks.gt.Ns) ks = Ns
      if (ks.lt.1) ks = 1                           
      Rclust = dss+ran2(idum)*(Rclmax-dss)
      Rclust2 = Rclust*Rclust
      attcln = attcln+1.
c      write(*,*) 'attclp = ',attclp
      call cln
      if (idonk.eq.1) goto 727
      dut = deltaur+deltauk
      if (dut.gt.100.d0) goto 727      
      if (dut.lt.0.d0) goto 7165
      if (dexp(-dut).lt.ran2(idum)) goto 727
 7165 acccln = acccln+1.      
      ur = ur+deltaur
      xs(ks) = txs
      ys(ks) = tys
      zs(ks) = tzs
      do kk = 1,kcln
      km = iclustn(kk)
      xs(km) = xcn(km)
      ys(km) = ycn(km)
      zs(km) = zcn(km)
      enddo
      do kk = 1,kclp
      km = iclustp(kk)
      xm(km) = xcp(km)
      ym(km) = ycp(km)
      zm(km) = zcp(km)
      enddo
      uk = uk+deltauk
      do kntot = 1,kntmax
      vcsum(kntot) = tvcsum(kntot)
      vssum(kntot) = tvssum(kntot)
      enddo      
      goto 6
 727  rejcln = rejcln+1.
      goto 6
      
      else
      atts = atts+1.
      call strans
      if (idonk.eq.1) goto 703
c      deltauk = 0.5d0*deltauk
c      deltaur = 0.5d0*deltaur
c      deltaurlj = urljn-urljo      
      deltau = deltauk+deltaur
      if (dexp(-deltau).lt.ran2(idum)) goto 703
      accs = accs+1.
      uk = uk+deltauk
      ur = ur+deltaur
c      urlj = urlj+deltaurlj
      xs(ip) = tx
      ys(ip) = ty
      zs(ip) = tz
      do kntot = 1,kntmax
      vcsum(kntot) = tvcsum(kntot)
      vssum(kntot) = tvssum(kntot)
      enddo
      goto 705
 703  rejs = rejs+1.
 705  goto 6
      endif
      endif
      endif
 6    continue
      if (mod(ksub,knsamp).eq.0) then 
      asampnz = asampnz+1.d0
c      avsdens = sdens+avsdens   
      avNs = rNs+avNs
      avNm = rNm+avNm  
      if (Ns.gt.maxNn) maxNn = Ns
      if (Nm.gt.maxNp) maxNp = Nm
      if (Ns.lt.minNn) minNn = Ns
      if (Nm.lt.minNp) minNp = Nm            
      rNmdist(Nm) = rNmdist(Nm)+1.
      rNsdist(Ns) = rNsdist(Ns)+1.      

      do i = 1,Nm
      tx = xm(i)
      ty = ym(i)
      tz = zm(i)
      idw = int((tz+ed)*rdeltaz)+1      
      rhozm(idw) = rhozm(idw)+1.
      do k = i+1,Nm
      rr = dist(tx,ty,tz,xm(k),ym(k),zm(k))
      if (rr.lt.ed) then
      idr = int(rr*rdr)+1
      grcc(idr) = grcc(idr)+1.
      endif
      enddo
      do k = 1,Ns
      rr = dist(tx,ty,tz,xs(k),ys(k),zs(k))
      if (rr.lt.ed) then
      idr = int(rr*rdr)+1
      grca(idr) = grca(idr)+1.
      endif
      enddo            
      enddo
      
      do i = 1,Ns
      tx = xs(i)
      ty = ys(i)
      tz = zs(i)
      idw = int((tz+ed)*rdeltaz)+1      
      rhozs(idw) = rhozs(idw)+1.
      do k = i+1,Ns
      rr = dist(tx,ty,tz,xs(k),ys(k),zs(k))
      if (rr.lt.ed) then
      idr = int(rr*rdr)+1
      graa(idr) = graa(idr)+1.
      endif
      enddo
      enddo

      endif
 276  continue
      suk = uk
      sur = ur
      suself = uself
c      surlj = urlj
      write(*,*)
      write(*,*) 'kmac = ',kmac
      write(*,*) 'suk,sur = ',suk,sur
      write(*,*) 'suself = ',suself
      call kvecs
      call check
      sfactor = (rNm*rnmval**2+rNs*rnsval**2)      
      uself = -sfactor*rbeta*alpha/dsqrt(pi)
c      uself = 0.5d0*uself         
      write(*,*) 'uk,ur = ',uk,ur
      write(*,*) 'uself = ',uself
      duk = dabs((suk-uk)/suk)
      dur = dabs((sur-ur)/sur)
      duself = dabs((suself-uself)/suself)
      write(*,*) 'duk = ',duk
      write(*,*) 'dur = ',dur
      write(*,*) 'duself = ',duself
      write(*,*) 'accm/attm = ',accm/attm
      write(*,*) 'accs/atts = ',accs/atts
      if (attcrem.gt.0.1) write(*,*)'acccrem/attcrem = ',acccrem/attcrem      
      if (attdesm.gt.0.1) write(*,*)'accdesm/attdesm = ',accdesm/attdesm
      if (attcres.gt.0.1) write(*,*)'acccres/attcres = ',acccres/attcres
      if (attdess.gt.0.1) write(*,*)'accdess/attdess = ',accdess/attdess             
      write(*,*) 'uk+ur = ',uk+ur 
      write(*,*) 'utot = ',uk+ur+uself
      write(*,*)
      write(*,*) 'faccclp = ',accclp/attclp
      write(*,*) 'facccln = ',acccln/attcln
      write(*,*) attclp,accclp,rejclp
      write(*,*) attcln,acccln,rejcln
      write(*,*)      
      write(*,*) 'rNs,avNs = ',rNs,avNs/asampnz
      write(*,*) 'rNm,avNm = ',rNm,avNm/asampnz
      cavsdens = 0.5d0*red2**2*(avNs-rnmval*avNm)/asampnz      
      write(*,*) 'avsdens = ',avsdens/sampsdens,cavsdens
      write(*,*)
c      call flush
 5    continue
      write(*,*) 
      rewind lll
      do i = 1,Nm
      write(lll,'(3f25.14)') xm(i),ym(i),zm(i)
      enddo
      rewind kkk
      do i = 1,Ns
      write(kkk,'(3f25.14)') xs(i),ys(i),zs(i)
      enddo
      write(*,*) 'asampnz = ',asampnz

      ih = nint(h/deltaz)
      write(*,*) 'asampnz,dV = ',asampnz,dV
      avnp = 0.d0
      avnn = 0.d0
      rewind inzp
      rewind inzn
      zz = h-halfh-0.5d0*deltaz
      do i = 1,ih
      zz = zz+deltaz 
      densm = rhozm(i)/(asampnz*dV)
      denss = rhozs(i)/(asampnz*dV)
      avnp = densm+avnp
      avnn = denss+avnn
      write(inzp,*) zz,densm
      write(inzn,*) zz,denss      
      enddo      

      rewind islump
      write(islump,*) idum
      rewind inf
      write(inf,*) Ns
      write(inf,*) Nm,nmval
      write(inf,*) chemp
      write(inf,*) donnan

      rewind Nnd
      rewind Npd
      write(*,*) 'minNn,maxNn = ',minNn,maxNn
      write(*,*) 'minNp,maxNp = ',minNp,maxNp
      do i = minNn,maxNn
      write(Nnd,*) i,rNsdist(i)/asampnz
      enddo
      do i = minNp,maxNp
      write(Npd,*) i,rNmdist(i)/asampnz
      enddo           

      dV = ed2*ed2*ed2
      write(*,*) 'asampnz = ',asampnz,deltaz
      nr = int(ed/deltaz)
      avNm = avNm/asampnz
      avNs = avNs/asampnz
      nr = nint(ed/deltaz)
      refmm = 0.5d0*avNm*(avNm-1.)*asampnz/ed2**3
      refms = avNm*avNs*asampnz/ed2**3
      refss = 0.5d0*avNs*(avNs-1.)*asampnz/ed2**3
      write(*,*) grcc(1),grcc(2),grcc(3)
      r = -deltaz*0.5d0
      rewind irdf
      do j = 1,nr
      r = r+deltaz
      rup = r+0.5d0*deltaz
      rlow = r-0.5d0*deltaz
      dV = 4.d0*pi*(rup**3-rlow**3)/3.d0
      refnomm = refmm*dV
      refnoms = refms*dV
      refnoss = refss*dV      
      rdfmm = grcc(j)/refnomm
      rdfms = grca(j)/refnoms
      rdfss = graa(j)/refnoss
      write(irdf,'(5e16.7)') r,rdfmm,rdfms,rdfss
      enddo

      ia =  70
      open (ia,file='a.pqr',form='formatted')
      rewind ia
      do i = 1,Nm
      write(ia,'(A,i7,A,i6,1f12.2,2f9.2,2f6.2)')
     *'ATOM',i,' unk unk ',i,zm(i),ym(i),xm(i),3.d0,2.d0
      enddo
      do i = 1,Ns
      write(ia,'(A,i7,A,i6,1f12.2,2f9.2,2f6.2)')
     *'ATOM',i,' uns uns ',i,zs(i),ys(i),xs(i),-1.d0,2.d0
      enddo

      stop
      end

      subroutine ksetup
      implicit double precision (a-h,o-z)   
      include 'yy.inc'                       
      zsum = 0.d0
      kntot = 0
      nx = 0
      xrk = dfloat(nx)*rkfact
      do ny = -kmax,kmax
      yrk = dfloat(ny)*rkfact         
      rhosq = xrk**2+yrk**2
      do nz = -kzmax,kzmax
      nsq = nx*nx+ny*ny+nz*nz         
      if (nsq.eq.0) goto 123

      if (nsq.gt.ksqmax) goto 123
      
      kntot = kntot+1
      zrk = dfloat(nz)*rkzfact               
      vxrk(kntot) = xrk
      vyrk(kntot) = yrk
      vzrk(kntot) = zrk
      rksq = rhosq+zrk**2
      fvec(kntot) = dexp(-rksq*efact)/rksq
 123  continue
      enddo
      enddo
      knt0 = kntot
      do nx = 1,kmax      
      xrk = dfloat(nx)*rkfact
      do ny = -kmax,kmax
      yrk = dfloat(ny)*rkfact         
      rhosq = xrk**2+yrk**2
      do nz = -kzmax,kzmax

      nsq = nx*nx+ny*ny+nz*nz
      if (nsq.gt.ksqmax) goto 125      
      
      kntot = kntot+1
      zrk = dfloat(nz)*rkzfact               
      vxrk(kntot) = xrk
      vyrk(kntot) = yrk
      vzrk(kntot) = zrk
      rksq = rhosq+zrk**2
      fvec(kntot) = dexp(-rksq*efact)/rksq
 125  continue      
      enddo
      enddo
      enddo
      kntmax = kntot
      write(*,*) 'kntmax = ',kntmax
      return
      end


      subroutine kvecs
      implicit double precision (a-h,o-z)   
      include 'yy.inc'                             
      zsum = 0.d0
      do kntot = 1,knt0
      xrk = vxrk(kntot)   
      yrk = vyrk(kntot)   
      zrk = vzrk(kntot)   
      csum = 0.d0
      ssum = 0.d0
      do i = 1,Nm
      csum = csum+dcos(xrk*xm(i)+yrk*ym(i)+zrk*zm(i))
      ssum = ssum+dsin(xrk*xm(i)+yrk*ym(i)+zrk*zm(i))
c      csum = csum-dcos(xrk*xm(i)+yrk*ym(i)+zrk*rzm(i))
c      ssum = ssum-dsin(xrk*xm(i)+yrk*ym(i)+zrk*rzm(i))      
      enddo
      csum = csum*rnmval
      ssum = ssum*rnmval
      do i = 1,Ns
      csum = csum-dcos(xrk*xs(i)+yrk*ys(i)+zrk*zs(i))
      ssum = ssum-dsin(xrk*xs(i)+yrk*ys(i)+zrk*zs(i))
c      csum = csum+dcos(xrk*xs(i)+yrk*ys(i)+zrk*rzs(i))
c      ssum = ssum+dsin(xrk*xs(i)+yrk*ys(i)+zrk*rzs(i))      
      enddo
      trhosq = csum**2+ssum**2
      zsum = zsum+fvec(kntot)*trhosq
      vcsum(kntot) = csum
      vssum(kntot) = ssum
      enddo
      sum = 0.d0
      do kntot = knt0+1,kntmax
      xrk = vxrk(kntot)   
      yrk = vyrk(kntot)   
      zrk = vzrk(kntot)   
      csum = 0.d0
      ssum = 0.d0
      do i = 1,Nm
      csum = csum+dcos(xrk*xm(i)+yrk*ym(i)+zrk*zm(i))
      ssum = ssum+dsin(xrk*xm(i)+yrk*ym(i)+zrk*zm(i))
c      csum = csum-dcos(xrk*xm(i)+yrk*ym(i)+zrk*rzm(i))
c      ssum = ssum-dsin(xrk*xm(i)+yrk*ym(i)+zrk*rzm(i))      
      enddo
      csum = csum*rnmval
      ssum = ssum*rnmval
      do i = 1,Ns
      csum = csum-dcos(xrk*xs(i)+yrk*ys(i)+zrk*zs(i))
      ssum = ssum-dsin(xrk*xs(i)+yrk*ys(i)+zrk*zs(i))
c      csum = csum+dcos(xrk*xs(i)+yrk*ys(i)+zrk*rzs(i))
c      ssum = ssum+dsin(xrk*xs(i)+yrk*ys(i)+zrk*rzs(i))      
      enddo
      trhosq = csum**2+ssum**2
      sum = sum+fvec(kntot)*trhosq
      vcsum(kntot) = csum
      vssum(kntot) = ssum
      enddo
      sum = 2.d0*sum+zsum
      sum = twopi*sum*rvolume
      uk = sum*bjerrum
      return
      end

      subroutine cln
      implicit double precision (a-h,o-z)      
      include 'yy.inc'                              
      idonk = 0
      dpx = cldps*(ran2(idum)-0.5d0)
      dpy = cldps*(ran2(idum)-0.5d0)
      dpz = cldps*(ran2(idum)-0.5d0)
      txs = xs(ks)+dpx
      tys = ys(ks)+dpy
      tzs = zs(ks)+dpz
      if (txs.gt.ed) then
      txs = txs-ed2
      elseif (txs.lt.-ed) then
      txs = txs+ed2
      endif
      if (tys.gt.ed) then
      tys = tys-ed2
      elseif (tys.lt.-ed) then
      tys = tys+ed2
      endif
      if (tzs.gt.ed) then
      tzs = tzs-ed2
      elseif (tzs.lt.-ed) then
      tzs = tzs+ed2
      endif
      kclp = 0
      ttxs = xs(ks)
      ttys = ys(ks)
      ttzs = zs(ks)
      udo = 0.d0
      do k = 1,Nm
      xcp(k) = xm(k)
      ycp(k) = ym(k)
      zcp(k) = zm(k)
      dx = dabs(ttxs-xm(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(ttys-ym(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(ttzs-zm(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.Rclust2) then
      kclp = kclp+1
      iclustp(kclp) = k         
      xx = xm(k)+dpx
      if (xx.gt.ed) then
      xx = xx-ed2
      elseif (xx.lt.-ed) then
      xx = xx+ed2
      endif
      xcp(k) = xx
      yy = ym(k)+dpy
      if (yy.gt.ed) then
      yy = yy-ed2
      elseif (yy.lt.-ed) then
      yy = yy+ed2
      endif
      ycp(k) = yy
      zz = zm(k)+dpz
      if (zz.gt.ed) then
      zz = zz-ed2
      elseif (zz.lt.-ed) then
      zz = zz+ed2
      endif
      zcp(k) = zz
      else
c     udo = -rnval/dsqrt(rmm2)+udo

      if (rmm2.lt.edsq) then
      tr = dsqrt(rmm2)
      udo = -rnval*terfc(alpha*tr)/tr+udo
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udo = -rnval*(dielc/tdielc-1.d0)/tr+udo
      endif
      endif 
      
      endif
      enddo
      
      kcln = 0
      do k = 1,ks-1
      xcn(k) = xs(k)
      ycn(k) = ys(k)
      zcn(k) = zs(k)
      dx = dabs(ttxs-xs(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(ttys-ys(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(ttzs-zs(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.Rclust2) then
      kcln = kcln+1
      iclustn(kcln) = k         
      xx = xs(k)+dpx
      if (xx.gt.ed) then
      xx = xx-ed2
      elseif (xx.lt.-ed) then
      xx = xx+ed2
      endif
      xcn(k) = xx
      yy = ys(k)+dpy
      if (yy.gt.ed) then
      yy = yy-ed2
      elseif (yy.lt.-ed) then
      yy = yy+ed2
      endif
      ycn(k) = yy
      zz = zs(k)+dpz
      if (zz.gt.ed) then
      zz = zz-ed2
      elseif (zz.lt.-ed) then
      zz = zz+ed2
      endif
      zcn(k) = zz

      else
c     udo = rnval2/dsqrt(rmm2)+udo

      if (rmm2.lt.edsq) then         
      tr = dsqrt(rmm2)
      udo = terfc(alpha*tr)/tr+udo
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udo = (dielc/tdielc-1.d0)/tr+udo
      endif      
      endif
      
      endif
      enddo

      do k = ks+1,Ns
      xcn(k) = xs(k)
      ycn(k) = ys(k)
      zcn(k) = zs(k)
      dx = dabs(ttxs-xs(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(ttys-ys(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(ttzs-zs(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.Rclust2) then
      kcln = kcln+1
      iclustn(kcln) = k         
      xx = xs(k)+dpx
      if (xx.gt.ed) then
      xx = xx-ed2
      elseif (xx.lt.-ed) then
      xx = xx+ed2
      endif
      xcn(k) = xx
      yy = ys(k)+dpy
      if (yy.gt.ed) then
      yy = yy-ed2
      elseif (yy.lt.-ed) then
      yy = yy+ed2
      endif
      ycn(k) = yy
      zz = zs(k)+dpz
      if (zz.gt.ed) then
      zz = zz-ed2
      elseif (zz.lt.-ed) then
      zz = zz+ed2
      endif
      zcn(k) = zz

      else
c     udo = rnval2/dsqrt(rmm2)+udo

      if (rmm2.lt.edsq) then         
      tr = dsqrt(rmm2)
      udo = terfc(alpha*tr)/tr+udo
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udo = (dielc/tdielc-1.d0)/tr+udo
      endif      
      endif
      
      endif
      enddo

      udn = 0.d0
      kkclp = 0
      do k = 1,Nm
      dx = dabs(txs-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(tys-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(tzs-zcp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dms2) goto 927
      if (rmm2.lt.Rclust2) then
      kkclp = kkclp+1
      else
c     udn = -rnval/dsqrt(rmm2)+udn

      if (rmm2.lt.edsq) then         
      tr = dsqrt(rmm2)
      udn = -rnval*terfc(alpha*tr)/tr+udn
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udn = -rnval*(dielc/tdielc-1.d0)/tr+udn
      endif      
      endif
      
      endif   
      enddo
      
      kkcln = 0
      do k = 1,ks-1
      dx = dabs(txs-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(tys-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(tzs-zcn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dmm2) goto 927
      if (rmm2.lt.Rclust2) then
      kkcln = kkcln+1
      else
c     udn = rnval2/dsqrt(rmm2)+udn

      if (rmm2.lt.edsq) then         
      tr = dsqrt(rmm2)
      udn = terfc(alpha*tr)/tr+udn
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udn = (dielc/tdielc-1.d0)/tr+udn
      endif      
      endif
      
      endif   
      enddo
      do k = ks+1,Ns
      dx = dabs(txs-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(tys-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(tzs-zcn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dmm2) goto 927
      if (rmm2.lt.Rclust2) then
      kkcln = kkcln+1
      else
c     udn = rnval2/dsqrt(rmm2)+udn

      if (rmm2.lt.edsq) then         
      tr = dsqrt(rmm2)
      udn = terfc(alpha*tr)/tr+udn
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udn = (dielc/tdielc-1.d0)/tr+udn
      endif      
      endif
      
      endif   
      enddo
      
      if (kkcln.ne.kcln.or.kkclp.ne.kclp) goto 927      
c     no cluster overlaps => calculate 
c     present and trial energies with other cluster members
      xcn(ks) = txs
      ycn(ks) = tys
      zcn(ks) = tzs      
      do kk = 1,kcln
      km = iclustn(kk)
      stx = xcn(km)
      sty = ycn(km)
      stz = zcn(km)
      ttx = xs(km)
      tty = ys(km)
      ttz = zs(km)
      
      do k = 1,km-1
      dx = abs(stx-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dss2) goto 927
c     udn = 1.d0/dsqrt(rmm2)+udn

      if (rmm2.lt.edsq) then
      tr = dsqrt(rmm2)
      udn = terfc(alpha*tr)/tr+udn
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udn = (dielc/tdielc-1.d0)/tr+udn
      endif      
      endif
      
      udo = uc(ttx,tty,ttz,xs(k),ys(k),zs(k))+udo      
      enddo
      do k = km+1,Ns
      dx = abs(stx-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dss2) goto 927
c     udn = 1.d0/dsqrt(rmm2)+udn

      if (rmm2.lt.edsq) then
      tr = dsqrt(rmm2)
      udn = terfc(alpha*tr)/tr+udn
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udn = (dielc/tdielc-1.d0)/tr+udn
      endif      
      endif

      udo = uc(ttx,tty,ttz,xs(k),ys(k),zs(k))+udo            
      enddo

      do k = 1,Nm      
      dx = abs(stx-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dms2) goto 927
c     udn = -rnval/dsqrt(rmm2)+udn

      if (rmm2.lt.edsq) then      
      tr = dsqrt(rmm2)
      udn = -rnval*terfc(alpha*tr)/tr+udn
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udn = -rnval*(dielc/tdielc-1.d0)/tr+udn
      endif            
      endif
      
      udo = -rnval*uc(ttx,tty,ttz,xm(k),ym(k),zm(k))+udo       
      enddo
      enddo

      do kk = 1,kclp
      km = iclustp(kk)
      stx = xcp(km)
      sty = ycp(km)
      stz = zcp(km)
      ttx = xm(km)
      tty = ym(km)
      ttz = zm(km)      
      do k = 1,km-1
      dx = abs(stx-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dmm2) goto 927
c     udn = rnval2/dsqrt(rmm2)+udn

      if (rmm2.lt.edsq) then      
      tr = dsqrt(rmm2)
      udn = rnval2*terfc(alpha*tr)/tr+udn
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udn = rnval2*(dielc/tdielc-1.d0)/tr+udn
      endif                  
      endif
      
      udo = rnval2*uc(ttx,tty,ttz,xm(k),ym(k),zm(k))+udo      
      enddo
      do k = km+1,Nm
      dx = abs(stx-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dmm2) goto 927
c     udn = rnval2/dsqrt(rmm2)+udn

      if (rmm2.lt.edsq) then      
      tr = dsqrt(rmm2)
      udn = rnval2*terfc(alpha*tr)/tr+udn
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udn = rnval2*(dielc/tdielc-1.d0)/tr+udn
      endif                  
      endif
      
      udo = rnval2*uc(ttx,tty,ttz,xm(k),ym(k),zm(k))+udo            
      enddo

      do k = 1,Ns
      dx = abs(stx-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dms2) goto 927
c     udn = -rnval/dsqrt(rmm2)+udn

      if (rmm2.lt.edsq) then      
      tr = dsqrt(rmm2)
      udn = -rnval*terfc(alpha*tr)/tr+udn
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udn = -rnval*(dielc/tdielc-1.d0)/tr+udn
      endif                  
      endif
      
      udo = -rnval*uc(ttx,tty,ttz,xs(k),ys(k),zs(k))+udo       
      enddo
      enddo
      deltaur = bjerrum*(udn-udo)      

      deltaukz = 0.d0      
      do kntot = 1,knt0
      xrk = vxrk(kntot)   
      yrk = vyrk(kntot)   
      zrk = vzrk(kntot)

      tx = xcn(ks)
      ty = ycn(ks)
      tz = zcn(ks)
      ttx = xs(ks)
      tty = ys(ks)
      ttz = zs(ks)                  
      targ = xrk*tx+yrk*ty+zrk*tz
      ttarg = xrk*ttx+yrk*tty+zrk*ttz
      dcsum = dcos(ttarg)-dcos(targ)
      dssum = dsin(ttarg)-dsin(targ)

      do kk = 1,kclp
      km = iclustp(kk)
      tx = xcp(km)
      ty = ycp(km)
      tz = zcp(km)
      ttx = xm(km)
      tty = ym(km)
      ttz = zm(km)                  
      targ = xrk*tx+yrk*ty+zrk*tz
      ttarg = xrk*ttx+yrk*tty+zrk*ttz
      dcsum = (-dcos(ttarg)+dcos(targ))*rnval+dcsum
      dssum = (-dsin(ttarg)+dsin(targ))*rnval+dssum
      enddo

      do kk = 1,kcln
      km = iclustn(kk)
      tx = xcn(km)
      ty = ycn(km)
      tz = zcn(km)
      ttx = xs(km)
      tty = ys(km)
      ttz = zs(km)
      targ = xrk*tx+yrk*ty+zrk*tz
      ttarg = xrk*ttx+yrk*tty+zrk*ttz
      dcsum = dcos(ttarg)-dcos(targ)+dcsum
      dssum = dsin(ttarg)-dsin(targ)+dssum
      enddo
      tvcsum(kntot) = vcsum(kntot)+dcsum
      tvssum(kntot) = vssum(kntot)+dssum      

      dc = tvcsum(kntot)**2-vcsum(kntot)**2
      ds = tvssum(kntot)**2-vssum(kntot)**2
      deltaukz = fvec(kntot)*(dc+ds)+deltaukz      
      enddo
            
      deltauk = 0.d0
      do kntot = knt0+1,kntmax
      xrk = vxrk(kntot)   
      yrk = vyrk(kntot)   
      zrk = vzrk(kntot)

      tx = xcn(ks)
      ty = ycn(ks)
      tz = zcn(ks)
      ttx = xs(ks)
      tty = ys(ks)
      ttz = zs(ks)                  
      targ = xrk*tx+yrk*ty+zrk*tz
      ttarg = xrk*ttx+yrk*tty+zrk*ttz
      dcsum = dcos(ttarg)-dcos(targ)
      dssum = dsin(ttarg)-dsin(targ)

      do kk = 1,kclp
      km = iclustp(kk)
      tx = xcp(km)
      ty = ycp(km)
      tz = zcp(km)
      ttx = xm(km)
      tty = ym(km)
      ttz = zm(km)                  
      targ = xrk*tx+yrk*ty+zrk*tz
      ttarg = xrk*ttx+yrk*tty+zrk*ttz
      dcsum = (-dcos(ttarg)+dcos(targ))*rnval+dcsum
      dssum = (-dsin(ttarg)+dsin(targ))*rnval+dssum
      enddo

      do kk = 1,kcln
      km = iclustn(kk)
      tx = xcn(km)
      ty = ycn(km)
      tz = zcn(km)
      ttx = xs(km)
      tty = ys(km)
      ttz = zs(km)
      targ = xrk*tx+yrk*ty+zrk*tz
      ttarg = xrk*ttx+yrk*tty+zrk*ttz
      dcsum = dcos(ttarg)-dcos(targ)+dcsum
      dssum = dsin(ttarg)-dsin(targ)+dssum
      enddo
      tvcsum(kntot) = vcsum(kntot)+dcsum
      tvssum(kntot) = vssum(kntot)+dssum      

      dc = tvcsum(kntot)**2-vcsum(kntot)**2
      ds = tvssum(kntot)**2-vssum(kntot)**2
      deltauk = fvec(kntot)*(dc+ds)+deltauk     
      enddo
      
      deltauk = 2.d0*deltauk+deltaukz
      deltauk = twopi*deltauk*rvolume
      deltauk = deltauk*bjerrum
      
      goto 929
 927  idonk = 1
 929  continue
      return
      end      
      
      subroutine clp
      implicit double precision (a-h,o-z)      
      include 'yy.inc'                              
      idonk = 0
      dpx = cldps*(ran2(idum)-0.5d0)
      dpy = cldps*(ran2(idum)-0.5d0)
      dpz = cldps*(ran2(idum)-0.5d0)
      txs = xm(ks)+dpx
      tys = ym(ks)+dpy
      tzs = zm(ks)+dpz
      if (txs.gt.ed) then
      txs = txs-ed2
      elseif (txs.lt.-ed) then
      txs = txs+ed2
      endif
      if (tys.gt.ed) then
      tys = tys-ed2
      elseif (tys.lt.-ed) then
      tys = tys+ed2
      endif
      if (tzs.gt.ed) then
      tzs = tzs-ed2
      elseif (tzs.lt.-ed) then
      tzs = tzs+ed2
      endif
      kcln = 0
      ttxs = xm(ks)
      ttys = ym(ks)
      ttzs = zm(ks)
      udo = 0.d0
      do k = 1,Ns
      xcn(k) = xs(k)
      ycn(k) = ys(k)
      zcn(k) = zs(k)
      dx = dabs(ttxs-xs(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(ttys-ys(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(ttzs-zs(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.Rclust2) then
      kcln = kcln+1
      iclustn(kcln) = k         
      xx = xs(k)+dpx
      if (xx.gt.ed) then
      xx = xx-ed2
      elseif (xx.lt.-ed) then
      xx = xx+ed2
      endif
      xcn(k) = xx
      yy = ys(k)+dpy
      if (yy.gt.ed) then
      yy = yy-ed2
      elseif (yy.lt.-ed) then
      yy = yy+ed2
      endif
      ycn(k) = yy
      zz = zs(k)+dpz
      if (zz.gt.ed) then
      zz = zz-ed2
      elseif (zz.lt.-ed) then
      zz = zz+ed2
      endif
      zcn(k) = zz
      else
c     udo = -rnval/dsqrt(rmm2)+udo

      if (rmm2.lt.edsq) then         
      tr = dsqrt(rmm2)
      udo = -rnval*terfc(alpha*tr)/tr+udo
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udo = -rnval*(dielc/tdielc-1.d0)/tr+udo
      endif                  
      endif
      
      endif
      enddo
      
      kclp = 0
      do k = 1,ks-1
      xcp(k) = xm(k)
      ycp(k) = ym(k)
      zcp(k) = zm(k)
      dx = dabs(ttxs-xm(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(ttys-ym(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(ttzs-zm(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.Rclust2) then
      kclp = kclp+1
      iclustp(kclp) = k         
      xx = xm(k)+dpx
      if (xx.gt.ed) then
      xx = xx-ed2
      elseif (xx.lt.-ed) then
      xx = xx+ed2
      endif
      xcp(k) = xx
      yy = ym(k)+dpy
      if (yy.gt.ed) then
      yy = yy-ed2
      elseif (yy.lt.-ed) then
      yy = yy+ed2
      endif
      ycp(k) = yy
      zz = zm(k)+dpz
      if (zz.gt.ed) then
      zz = zz-ed2
      elseif (zz.lt.-ed) then
      zz = zz+ed2
      endif
      zcp(k) = zz

      else
c     udo = rnval2/dsqrt(rmm2)+udo

      if (rmm2.lt.edsq) then         
      tr = dsqrt(rmm2)
      udo = rnval2*terfc(alpha*tr)/tr+udo
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udo = rnval2*(dielc/tdielc-1.d0)/tr+udo
      endif                        
      endif
      
      endif
      enddo

      do k = ks+1,Nm
      xcp(k) = xm(k)
      ycp(k) = ym(k)
      zcp(k) = zm(k)
      dx = dabs(ttxs-xm(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(ttys-ym(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(ttzs-zm(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.Rclust2) then
      kclp = kclp+1
      iclustp(kclp) = k         
      xx = xm(k)+dpx
      if (xx.gt.ed) then
      xx = xx-ed2
      elseif (xx.lt.-ed) then
      xx = xx+ed2
      endif
      xcp(k) = xx
      yy = ym(k)+dpy
      if (yy.gt.ed) then
      yy = yy-ed2
      elseif (yy.lt.-ed) then
      yy = yy+ed2
      endif
      ycp(k) = yy
      zz = zm(k)+dpz
      if (zz.gt.ed) then
      zz = zz-ed2
      elseif (zz.lt.-ed) then
      zz = zz+ed2
      endif
      zcp(k) = zz

      else
c     udo = rnval2/dsqrt(rmm2)+udo

      if (rmm2.lt.edsq) then         
      tr = dsqrt(rmm2)
      udo = rnval2*terfc(alpha*tr)/tr+udo
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udo = rnval2*(dielc/tdielc-1.d0)/tr+udo
      endif                        
      endif
      
      endif
      enddo            

      udn = 0.d0
      kkcln = 0
      do k = 1,Ns
      dx = dabs(txs-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(tys-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(tzs-zcn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dms2) goto 927
      if (rmm2.lt.Rclust2) then
      kkcln = kkcln+1
      else
c     udn = -rnval/dsqrt(rmm2)+udn

      if (rmm2.lt.edsq) then         
      tr = dsqrt(rmm2)
      udn = -rnval*terfc(alpha*tr)/tr+udn
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udn = -rnval*(dielc/tdielc-1.d0)/tr+udn
      endif                              
      endif
      
      endif   
      enddo

      kkclp = 0
      do k = 1,ks-1
      dx = dabs(txs-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(tys-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(tzs-zcp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dmm2) goto 927
      if (rmm2.lt.Rclust2) then
      kkclp = kkclp+1
      else
c     udn = rnval2/dsqrt(rmm2)+udn

      if (rmm2.lt.edsq) then         
      tr = dsqrt(rmm2)
      udn = rnval2*terfc(alpha*tr)/tr+udn
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udn = rnval2*(dielc/tdielc-1.d0)/tr+udn
      endif                                    
      endif      
      
      endif   
      enddo
      do k = ks+1,Nm
      dx = dabs(txs-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = dabs(tys-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = dabs(tzs-zcp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dmm2) goto 927
      if (rmm2.lt.Rclust2) then
      kkclp = kkclp+1
      else
c     udn = rnval2/dsqrt(rmm2)+udn

      if (rmm2.lt.edsq) then         
      tr = dsqrt(rmm2)
      udn = rnval2*terfc(alpha*tr)/tr+udn
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udn = rnval2*(dielc/tdielc-1.d0)/tr+udn
      endif                                    
      endif
      
      endif   
      enddo
      
      if (kkcln.ne.kcln.or.kkclp.ne.kclp) goto 927      
c     no cluster overlaps => calculate 
c     present and trial energies with other cluster members
      xcp(ks) = txs
      ycp(ks) = tys
      zcp(ks) = tzs      
      do kk = 1,kcln
      km = iclustn(kk)
      stx = xcn(km)
      sty = ycn(km)
      stz = zcn(km)
      ttx = xs(km)
      tty = ys(km)
      ttz = zs(km)
      
      do k = 1,km-1
      dx = abs(stx-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dss2) goto 927
c     udn = 1.d0/dsqrt(rmm2)+udn

      if (rmm2.lt.edsq) then      
      tr = dsqrt(rmm2)
      udn = terfc(alpha*tr)/tr+udn
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udn = (dielc/tdielc-1.d0)/tr+udn
      endif                                    
      endif
      
      udo = uc(ttx,tty,ttz,xs(k),ys(k),zs(k))+udo      
      enddo
      do k = km+1,Ns
      dx = abs(stx-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dss2) goto 927
c     udn = 1.d0/dsqrt(rmm2)+udn

      if (rmm2.lt.edsq) then      
      tr = dsqrt(rmm2)
      udn = terfc(alpha*tr)/tr+udn
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udn = (dielc/tdielc-1.d0)/tr+udn
      endif                                    
      endif
      
      udo = uc(ttx,tty,ttz,xs(k),ys(k),zs(k))+udo            
      enddo

      do k = 1,Nm      
      dx = abs(stx-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dms2) goto 927
c     udn = -rnval/dsqrt(rmm2)+udn

      if (rmm2.lt.edsq) then      
      tr = dsqrt(rmm2)
      udn = -rnval*terfc(alpha*tr)/tr+udn
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udn = -rnval*(dielc/tdielc-1.d0)/tr+udn
      endif                                    
      endif
      
      udo = -rnval*uc(ttx,tty,ttz,xm(k),ym(k),zm(k))+udo       
      enddo
      enddo

      do kk = 1,kclp
      km = iclustp(kk)
      stx = xcp(km)
      sty = ycp(km)
      stz = zcp(km)
      ttx = xm(km)
      tty = ym(km)
      ttz = zm(km)      
      do k = 1,km-1
      dx = abs(stx-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dmm2) goto 927
c     udn = rnval2/dsqrt(rmm2)+udn

      if (rmm2.lt.edsq) then      
      tr = dsqrt(rmm2)
      udn = rnval2*terfc(alpha*tr)/tr+udn
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udn = rnval2*(dielc/tdielc-1.d0)/tr+udn
      endif                                    
      endif
      
      udo = rnval2*uc(ttx,tty,ttz,xm(k),ym(k),zm(k))+udo      
      enddo
      do k = km+1,Nm
      dx = abs(stx-xcp(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycp(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcp(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dmm2) goto 927
c     udn = rnval2/dsqrt(rmm2)+udn

      if (rmm2.lt.edsq) then      
      tr = dsqrt(rmm2)
      udn = rnval2*terfc(alpha*tr)/tr+udn
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udn = rnval2*(dielc/tdielc-1.d0)/tr+udn
      endif                                    
      endif
      
      udo = rnval2*uc(ttx,tty,ttz,xm(k),ym(k),zm(k))+udo            
      enddo

      do k = 1,Ns
      dx = abs(stx-xcn(k))
      dx = dx-aint(dx*red)*ed2
      dy = abs(sty-ycn(k))
      dy = dy-aint(dy*red)*ed2
      dz = abs(stz-zcn(k))
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz
      if (rmm2.lt.dms2) goto 927
c     udn = -rnval/dsqrt(rmm2)+udn

      if (rmm2.lt.edsq) then      
      tr = dsqrt(rmm2)
      udn = -rnval*terfc(alpha*tr)/tr+udn
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      udn = -rnval*(dielc/tdielc-1.d0)/tr+udn
      endif                                    
      endif
      
      udo = -rnval*uc(ttx,tty,ttz,xs(k),ys(k),zs(k))+udo       
      enddo
      enddo
      deltaur = bjerrum*(udn-udo)      

      deltaukz = 0.d0      
      do kntot = 1,knt0
      xrk = vxrk(kntot)   
      yrk = vyrk(kntot)   
      zrk = vzrk(kntot)

      tx = xcp(ks)
      ty = ycp(ks)
      tz = zcp(ks)
      ttx = xm(ks)
      tty = ym(ks)
      ttz = zm(ks)                  
      targ = xrk*tx+yrk*ty+zrk*tz
      ttarg = xrk*ttx+yrk*tty+zrk*ttz
      dcsum = (-dcos(ttarg)+dcos(targ))*rnval
      dssum = (-dsin(ttarg)+dsin(targ))*rnval      

      do kk = 1,kclp
      km = iclustp(kk)
      tx = xcp(km)
      ty = ycp(km)
      tz = zcp(km)
      ttx = xm(km)
      tty = ym(km)
      ttz = zm(km)                  
      targ = xrk*tx+yrk*ty+zrk*tz
      ttarg = xrk*ttx+yrk*tty+zrk*ttz
      dcsum = (-dcos(ttarg)+dcos(targ))*rnval+dcsum
      dssum = (-dsin(ttarg)+dsin(targ))*rnval+dssum
      enddo

      do kk = 1,kcln
      km = iclustn(kk)
      tx = xcn(km)
      ty = ycn(km)
      tz = zcn(km)
      ttx = xs(km)
      tty = ys(km)
      ttz = zs(km)
      targ = xrk*tx+yrk*ty+zrk*tz
      ttarg = xrk*ttx+yrk*tty+zrk*ttz
      dcsum = dcos(ttarg)-dcos(targ)+dcsum
      dssum = dsin(ttarg)-dsin(targ)+dssum
      enddo
      tvcsum(kntot) = vcsum(kntot)+dcsum
      tvssum(kntot) = vssum(kntot)+dssum      

      dc = tvcsum(kntot)**2-vcsum(kntot)**2
      ds = tvssum(kntot)**2-vssum(kntot)**2
      deltaukz = fvec(kntot)*(dc+ds)+deltaukz      
      enddo
      
      
      deltauk = 0.d0
      do kntot = knt0+1,kntmax
      xrk = vxrk(kntot)   
      yrk = vyrk(kntot)   
      zrk = vzrk(kntot)

      tx = xcp(ks)
      ty = ycp(ks)
      tz = zcp(ks)
      ttx = xm(ks)
      tty = ym(ks)
      ttz = zm(ks)                  
      targ = xrk*tx+yrk*ty+zrk*tz
      ttarg = xrk*ttx+yrk*tty+zrk*ttz
      dcsum = (-dcos(ttarg)+dcos(targ))*rnval
      dssum = (-dsin(ttarg)+dsin(targ))*rnval      

      do kk = 1,kclp
      km = iclustp(kk)
      tx = xcp(km)
      ty = ycp(km)
      tz = zcp(km)
      ttx = xm(km)
      tty = ym(km)
      ttz = zm(km)                  
      targ = xrk*tx+yrk*ty+zrk*tz
      ttarg = xrk*ttx+yrk*tty+zrk*ttz
      dcsum = (-dcos(ttarg)+dcos(targ))*rnval+dcsum
      dssum = (-dsin(ttarg)+dsin(targ))*rnval+dssum
      enddo

      do kk = 1,kcln
      km = iclustn(kk)
      tx = xcn(km)
      ty = ycn(km)
      tz = zcn(km)
      ttx = xs(km)
      tty = ys(km)
      ttz = zs(km)
      targ = xrk*tx+yrk*ty+zrk*tz
      ttarg = xrk*ttx+yrk*tty+zrk*ttz
      dcsum = dcos(ttarg)-dcos(targ)+dcsum
      dssum = dsin(ttarg)-dsin(targ)+dssum
      enddo
      tvcsum(kntot) = vcsum(kntot)+dcsum
      tvssum(kntot) = vssum(kntot)+dssum      

      dc = tvcsum(kntot)**2-vcsum(kntot)**2
      ds = tvssum(kntot)**2-vssum(kntot)**2
      deltauk = fvec(kntot)*(dc+ds)+deltauk     
      enddo
      
      deltauk = 2.d0*deltauk+deltaukz
      deltauk = twopi*deltauk*rvolume
      deltauk = deltauk*bjerrum
      
      goto 929
 927  idonk = 1
 929  continue
      return
      end      
     
      subroutine crem
      implicit double precision (a-h,o-z)      
      include 'yy.inc'      
      tx = (ran2(idum)-0.5d0)*ed2
      ty = (ran2(idum)-0.5d0)*ed2
      tz = (ran2(idum)-0.5d0)*h
      
      urljn = 0.d0
      urn = 0.d0
      do j = 1,Nm

c     mA-mA:
      trho2 = rho2(tx,ty,xm(j),ym(j))
      dz = zm(j)-tz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2   
      if (tr2.lt.dmm2) goto 712
c     urljn = (dmm2/tr2)**6+urljn

      if (tr2.lt.edsq) then      
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urn = urn+terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urn = (dielc/tdielc-1.d0)*rtr+urn
      endif                                    
      endif
      
      enddo
      urn = urn*rnmval

      do j = 1,Ns
c     mA-sA:
      trho2 = rho2(tx,ty,xs(j),ys(j))
      dz = zs(j)-tz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2   
      if (tr2.lt.dms2) goto 712
c     urljn = (dms2/tr2)**6+urljn

      if (tr2.lt.edsq) then      
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urn = urn-terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urn = -(dielc/tdielc-1.d0)*rtr+urn
      endif                                          
      endif
      
      enddo
      
      urn = urn*rnmval
      deltaur = bjerrum*urn
      
      deltaukz = 0.d0
      do kntot = 1,knt0
      xrk = vxrk(kntot)   
      yrk = vyrk(kntot)   
      zrk = vzrk(kntot)
      targ = xrk*tx+yrk*ty+zrk*tz
c      rtarg = xrk*tx+yrk*ty+zrk*rtz
c      dcsum = dcos(targ)-dcos(rtarg)
c      dssum = dsin(targ)-dsin(rtarg)           
      dcsum = dcos(targ)
      dssum = dsin(targ)
      tvcsum(kntot) = vcsum(kntot)+dcsum*rnmval
      tvssum(kntot) = vssum(kntot)+dssum*rnmval      
      dc = tvcsum(kntot)**2-vcsum(kntot)**2
      ds = tvssum(kntot)**2-vssum(kntot)**2
      deltaukz = fvec(kntot)*(dc+ds)+deltaukz
      enddo

      deltauk = 0.d0
      do kntot = knt0+1,kntmax
      xrk = vxrk(kntot)   
      yrk = vyrk(kntot)   
      zrk = vzrk(kntot)
      targ = xrk*tx+yrk*ty+zrk*tz
c      rtarg = xrk*tx+yrk*ty+zrk*rtz
c      dcsum = dcos(targ)-dcos(rtarg)
c      dssum = dsin(targ)-dsin(rtarg)
      dcsum = dcos(targ)
      dssum = dsin(targ)
      tvcsum(kntot) = vcsum(kntot)+dcsum*rnmval
      tvssum(kntot) = vssum(kntot)+dssum*rnmval
      dc = tvcsum(kntot)**2-vcsum(kntot)**2
      ds = tvssum(kntot)**2-vssum(kntot)**2
      deltauk = fvec(kntot)*(dc+ds)+deltauk
      enddo
      deltauk = 2.d0*deltauk+deltaukz
      deltauk = twopi*deltauk*rvolume
      deltauk = deltauk*bjerrum
      
      goto 716
 712  idonk = 1
 716  continue      
      return
      end


      subroutine desm
      implicit double precision (a-h,o-z)      
      include 'yy.inc'      
      tx = xm(ks)
      ty = ym(ks)
      tz = zm(ks)

      urn = 0.d0
      urljn = 0.d0
      do j = 1,ks-1

c     mA-mA:
      trho2 = rho2(tx,ty,xm(j),ym(j))
      dz = zm(j)-tz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2   
c     if (tr2.lt.dmm2) goto 712
c      urljn = (dmm2/tr2)**6+urljn

      if (tr2.lt.edsq) then      
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urn = urn+terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urn = (dielc/tdielc-1.d0)*rtr+urn
      endif                                          
      endif
      
      enddo

      do j = ks+1,Nm

c     mA-mA:
      trho2 = rho2(tx,ty,xm(j),ym(j))
      dz = zm(j)-tz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2   
c     if (tr2.lt.dmm2) goto 712
c     urljn = (dmm2/tr2)**6+urljn

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urn = urn+terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urn = (dielc/tdielc-1.d0)*rtr+urn
      endif                                          
      endif
      
      enddo      

      urn = urn*rnmval

      do j = 1,Ns
c     mA-sA:
      trho2 = rho2(tx,ty,xs(j),ys(j))
      dz = zs(j)-tz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2   
c     if (tr2.lt.dms2) goto 712
c     urljn = (dms2/tr2)**6+urljn

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urn = urn-terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urn = -(dielc/tdielc-1.d0)*rtr+urn
      endif                                          
      endif
      
      enddo
      
      urn = urn*rnmval
      deltaur = bjerrum*urn
      
      deltaukz = 0.d0
      do kntot = 1,knt0
      xrk = vxrk(kntot)   
      yrk = vyrk(kntot)   
      zrk = vzrk(kntot)
      targ = xrk*tx+yrk*ty+zrk*tz
c      rtarg = xrk*tx+yrk*ty+zrk*rtz
c      dcsum = -dcos(targ)+dcos(rtarg)
c      dssum = -dsin(targ)+dsin(rtarg)
      dcsum = -dcos(targ)
      dssum = -dsin(targ)
      tvcsum(kntot) = vcsum(kntot)+dcsum*rnmval
      tvssum(kntot) = vssum(kntot)+dssum*rnmval
      dc = tvcsum(kntot)**2-vcsum(kntot)**2
      ds = tvssum(kntot)**2-vssum(kntot)**2
      deltaukz = fvec(kntot)*(dc+ds)+deltaukz
      enddo

      deltauk = 0.d0
      do kntot = knt0+1,kntmax
      xrk = vxrk(kntot)   
      yrk = vyrk(kntot)   
      zrk = vzrk(kntot)
      targ = xrk*tx+yrk*ty+zrk*tz
c      rtarg = xrk*tx+yrk*ty+zrk*rtz
c      dcsum = -dcos(targ)+dcos(rtarg)
c      dssum = -dsin(targ)+dsin(rtarg)
      dcsum = -dcos(targ)
      dssum = -dsin(targ)
      tvcsum(kntot) = vcsum(kntot)+dcsum*rnmval
      tvssum(kntot) = vssum(kntot)+dssum*rnmval
      dc = tvcsum(kntot)**2-vcsum(kntot)**2
      ds = tvssum(kntot)**2-vssum(kntot)**2
      deltauk = fvec(kntot)*(dc+ds)+deltauk
      enddo
      deltauk = 2.d0*deltauk+deltaukz
      deltauk = twopi*deltauk*rvolume
      deltauk = deltauk*bjerrum      
      return
      end


      subroutine cres
      implicit double precision (a-h,o-z)      
      include 'yy.inc'      
      tx = (ran2(idum)-0.5d0)*ed2
      ty = (ran2(idum)-0.5d0)*ed2
      tz = (ran2(idum)-0.5d0)*h      

      urn = 0.d0
      urljn = 0.d0      
      do j = 1,Nm
c     sA-mA         
      trho2 = rho2(tx,ty,xm(j),ym(j))
      dz = zm(j)-tz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2   
      if (tr2.lt.dms2) goto 712
c     urljn = (dms2/tr2)**6+urljn

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urn = urn-terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urn = -(dielc/tdielc-1.d0)*rtr+urn
      endif                                          
      endif
      
      enddo
      urn = urn*rnmval

      do j = 1,Ns
c     sA-sA
      trho2 = rho2(tx,ty,xs(j),ys(j))
      dz = zs(j)-tz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2   
      if (tr2.lt.dss2) goto 712
c     urljn = (dss2/tr2)**6+urljn

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urn = urn+terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urn = (dielc/tdielc-1.d0)*rtr+urn
      endif                                          
      endif
      
      enddo

      deltaur = bjerrum*urn

      deltaukz = 0.d0
      do kntot = 1,knt0
      xrk = vxrk(kntot)   
      yrk = vyrk(kntot)   
      zrk = vzrk(kntot)   
      targ = xrk*tx+yrk*ty+zrk*tz
c      rtarg = xrk*tx+yrk*ty+zrk*rtz
c      dcsum = -dcos(targ)+dcos(rtarg)
c      dssum = -dsin(targ)+dsin(rtarg)
      dcsum = -dcos(targ)
      dssum = -dsin(targ)
      tvcsum(kntot) = vcsum(kntot)+dcsum
      tvssum(kntot) = vssum(kntot)+dssum
      dc = tvcsum(kntot)**2-vcsum(kntot)**2
      ds = tvssum(kntot)**2-vssum(kntot)**2
      deltaukz = fvec(kntot)*(dc+ds)+deltaukz
      enddo

      deltauk = 0.d0
      do kntot = knt0+1,kntmax
      xrk = vxrk(kntot)   
      yrk = vyrk(kntot)   
      zrk = vzrk(kntot)
      targ = xrk*tx+yrk*ty+zrk*tz
c      rtarg = xrk*tx+yrk*ty+zrk*rtz
c      dcsum = -dcos(targ)+dcos(rtarg)
c      dssum = -dsin(targ)+dsin(rtarg)
      dcsum = -dcos(targ)
      dssum = -dsin(targ)
      tvcsum(kntot) = vcsum(kntot)+dcsum
      tvssum(kntot) = vssum(kntot)+dssum
      dc = tvcsum(kntot)**2-vcsum(kntot)**2
      ds = tvssum(kntot)**2-vssum(kntot)**2
      deltauk = fvec(kntot)*(dc+ds)+deltauk
      enddo
      deltauk = 2.d0*deltauk+deltaukz
      deltauk = twopi*deltauk*rvolume
      deltauk = deltauk*bjerrum

      goto 716
 712  idonk = 1
 716  continue
      return
      end

      subroutine dess
      implicit double precision (a-h,o-z)      
      include 'yy.inc'
      tx = xs(ks)
      ty = ys(ks)
      tz = zs(ks)

      urn = 0.d0
      urljn = 0.d0      
      do j = 1,Nm
c     sA-mA         
      trho2 = rho2(tx,ty,xm(j),ym(j))
      dz = zm(j)-tz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2   
c     if (tr2.lt.dms2) goto 712
c     urljn = (dms2/tr2)**6+urljn

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urn = urn-terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urn = -(dielc/tdielc-1.d0)*rtr+urn
      endif                                          
      endif
      
      enddo
      urn = urn*rnmval

      do j = 1,ks-1
c     sA-sA
      trho2 = rho2(tx,ty,xs(j),ys(j))
      dz = zs(j)-tz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2   
c     if (tr2.lt.dss2) goto 712
c     urljn = (dss2/tr2)**6+urljn

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urn = urn+terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urn = (dielc/tdielc-1.d0)*rtr+urn
      endif                                          
      endif
      
      enddo

      do j = ks+1,Ns
c     sA-sA
      trho2 = rho2(tx,ty,xs(j),ys(j))
      dz = zs(j)-tz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2   
c     if (tr2.lt.dss2) goto 712
c     urljn = (dss2/tr2)**6+urljn

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urn = urn+terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urn = (dielc/tdielc-1.d0)*rtr+urn
      endif                                          
      endif
      
      enddo      

      deltaur = bjerrum*urn

      deltaukz = 0.d0
      do kntot = 1,knt0
      xrk = vxrk(kntot)   
      yrk = vyrk(kntot)   
      zrk = vzrk(kntot)   
      targ = xrk*tx+yrk*ty+zrk*tz
c      rtarg = xrk*tx+yrk*ty+zrk*rtz
c      dcsum = dcos(targ)-dcos(rtarg)
c      dssum = dsin(targ)-dsin(rtarg)
      dcsum = dcos(targ)
      dssum = dsin(targ)
      tvcsum(kntot) = vcsum(kntot)+dcsum
      tvssum(kntot) = vssum(kntot)+dssum
      dc = tvcsum(kntot)**2-vcsum(kntot)**2
      ds = tvssum(kntot)**2-vssum(kntot)**2
      deltaukz = fvec(kntot)*(dc+ds)+deltaukz
      enddo

      deltauk = 0.d0
      do kntot = knt0+1,kntmax
      xrk = vxrk(kntot)   
      yrk = vyrk(kntot)   
      zrk = vzrk(kntot)
      targ = xrk*tx+yrk*ty+zrk*tz
c      rtarg = xrk*tx+yrk*ty+zrk*rtz
c      dcsum = dcos(targ)-dcos(rtarg)
c      dssum = dsin(targ)-dsin(rtarg)
      dcsum = dcos(targ)
      dssum = dsin(targ)
      tvcsum(kntot) = vcsum(kntot)+dcsum
      tvssum(kntot) = vssum(kntot)+dssum
      dc = tvcsum(kntot)**2-vcsum(kntot)**2
      ds = tvssum(kntot)**2-vssum(kntot)**2
      deltauk = fvec(kntot)*(dc+ds)+deltauk      
      enddo
      deltauk = 2.d0*deltauk+deltaukz
      deltauk = twopi*deltauk*rvolume
      deltauk = deltauk*bjerrum

      return
      end            
      

      subroutine mtrans
      implicit double precision (a-h,o-z)   
      include 'yy.inc'
      dumi = 0.d0
      idonk = 0
      ip =  int(1.d0+rNm*ran2(idum))
      if (ran2(idum).lt.fracslump) then
      tx = ed2*(ran2(idum)-0.5d0)
      ty = ed2*(ran2(idum)-0.5d0)                  
      tz = ed2*(ran2(idum)-0.5d0)
      else
      tx = xm(ip)+dpp*(ran2(idum)-0.5d0)
      ty = ym(ip)+dpp*(ran2(idum)-0.5d0)
      tx = tx-ed2*anint(tx*red2)
      ty = ty-ed2*anint(ty*red2)         
      tz = zm(ip)+dpp*(ran2(idum)-0.5d0)
      tz = tz-h*anint(tz*rh)               
      endif   
      ttx = xm(ip)
      tty = ym(ip)
      ttz = zm(ip)

      urold = 0.d0
      urnew = 0.d0
      umiold = 0.d0
      uminew = 0.d0
      urljn = 0.d0
      urljo = 0.d0 
      do j = 1,ip-1

c     mA-mA:
      trho2 = rho2(tx,ty,xm(j),ym(j))
      dz = zm(j)-tz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2   
      if (tr2.lt.dmm2) goto 712
c     urljn = (dmm2/tr2)**6+urljn

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urnew = urnew+terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urnew = (dielc/tdielc-1.d0)*rtr+urnew
      endif                                          
      endif
      
c      uminew = uminew+rtr

c     mA-mA:
      trho2 = rho2(ttx,tty,xm(j),ym(j))
      dz = zm(j)-ttz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2
c     urljo = (dmm2/tr2)**6+urljo

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urold = urold+terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urold = (dielc/tdielc-1.d0)*rtr+urold
      endif                                                
      endif
      
c      umiold = umiold+rtr
      enddo
      
      do j = ip+1,Nm
c     mA-mA:
      trho2 = rho2(tx,ty,xm(j),ym(j))
      dz = zm(j)-tz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2   
      if (tr2.lt.dmm2) goto 712
c     urljn = (dmm2/tr2)**6+urljn

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urnew = urnew+terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urnew = (dielc/tdielc-1.d0)*rtr+urnew
      endif                                                
      endif
      
c      uminew = uminew+rtr

c     mA-mA:
      trho2 = rho2(ttx,tty,xm(j),ym(j))
      dz = zm(j)-ttz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2
c     urljo = (dmm2/tr2)**6+urljo

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urold = urold+terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urold = (dielc/tdielc-1.d0)*rtr+urold
      endif                                                
      endif
      
c      umiold = umiold+rtr
      enddo
      urnew = urnew*rnmval
      urold = urold*rnmval
c      uminew = uminew*rnmval
c      umiold = umiold*rnmval      

      do j = 1,Ns
c     mA-sA:
      trho2 = rho2(tx,ty,xs(j),ys(j))
      dz = zs(j)-tz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2   
      if (tr2.lt.dms2) goto 712
c     urljn = (dms2/tr2)**6+urljn

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urnew = urnew-terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urnew = -(dielc/tdielc-1.d0)*rtr+urnew
      endif                                                
      endif
      
c      uminew = uminew-rtr

c     mA-sA:
      trho2 = rho2(ttx,tty,xs(j),ys(j))
      dz = zs(j)-ttz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2
c     urljo = (dms2/tr2)**6+urljo

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urold = urold-terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urold = -(dielc/tdielc-1.d0)*rtr+urold
      endif                                                
      endif
      
c      umiold = umiold-rtr
      enddo
      
      urnew = urnew*rnmval
      urold = urold*rnmval
c      uminew = uminew*rnmval
c      umiold = umiold*rnmval
c      dumi = bjerrum*(uminew-umiold)
      deltaur = bjerrum*(urnew-urold)
      
      deltaukz = 0.d0
      do kntot = 1,knt0
      xrk = vxrk(kntot)   
      yrk = vyrk(kntot)   
      zrk = vzrk(kntot)
      targ = xrk*tx+yrk*ty+zrk*tz
      ttarg = xrk*ttx+yrk*tty+zrk*ttz
c      rtarg = xrk*tx+yrk*ty+zrk*rtz
c      rttarg = xrk*ttx+yrk*tty+zrk*rttz      
c      dcsum = -dcos(ttarg)+dcos(targ)+dcos(rttarg)-dcos(rtarg)
c      dssum = -dsin(ttarg)+dsin(targ)+dsin(rttarg)-dsin(rtarg)
      dcsum = -dcos(ttarg)+dcos(targ)
      dssum = -dsin(ttarg)+dsin(targ)
      tvcsum(kntot) = vcsum(kntot)+dcsum*rnmval
      tvssum(kntot) = vssum(kntot)+dssum*rnmval
      dc = tvcsum(kntot)**2-vcsum(kntot)**2
      ds = tvssum(kntot)**2-vssum(kntot)**2
      deltaukz = fvec(kntot)*(dc+ds)+deltaukz
      enddo

      deltauk = 0.d0
      do kntot = knt0+1,kntmax
      xrk = vxrk(kntot)   
      yrk = vyrk(kntot)   
      zrk = vzrk(kntot)
      targ = xrk*tx+yrk*ty+zrk*tz
      ttarg = xrk*ttx+yrk*tty+zrk*ttz      
c      rtarg = xrk*tx+yrk*ty+zrk*rtz
c      rttarg = xrk*ttx+yrk*tty+zrk*rttz      
c      dcsum = -dcos(ttarg)+dcos(targ)+dcos(rttarg)-dcos(rtarg)
c      dssum = -dsin(ttarg)+dsin(targ)+dsin(rttarg)-dsin(rtarg)
      dcsum = -dcos(ttarg)+dcos(targ)
      dssum = -dsin(ttarg)+dsin(targ)
      tvcsum(kntot) = vcsum(kntot)+dcsum*rnmval
      tvssum(kntot) = vssum(kntot)+dssum*rnmval
      dc = tvcsum(kntot)**2-vcsum(kntot)**2
      ds = tvssum(kntot)**2-vssum(kntot)**2
      deltauk = fvec(kntot)*(dc+ds)+deltauk
      enddo
      deltauk = 2.d0*deltauk+deltaukz
      deltauk = twopi*deltauk*rvolume
      deltauk = deltauk*bjerrum
      goto 716
 712  idonk = 1
 716  continue
      return
      end

      subroutine strans
      implicit double precision (a-h,o-z)   
      include 'yy.inc'
      idonk = 0
      ip =  int(1.d0+rNs*ran2(idum))
      if (ran2(idum).lt.fracslump) then
      tx = ed2*(ran2(idum)-0.5d0)
      ty = ed2*(ran2(idum)-0.5d0)         
      tz = ed2*(ran2(idum)-0.5d0)
      else
      tx = xs(ip)+dps*(ran2(idum)-0.5d0)
      ty = ys(ip)+dps*(ran2(idum)-0.5d0)
      tx = tx-ed2*anint(tx*red2)
      ty = ty-ed2*anint(ty*red2)         
      tz = zs(ip)+dps*(ran2(idum)-0.5d0)
      tz = tz-h*anint(tz*rh)                
      endif   
      ttx = xs(ip)
      tty = ys(ip)
      ttz = zs(ip)
      urold = 0.d0
      urnew = 0.d0
      umiold = 0.d0
      uminew = 0.d0
      urljo = 0.d0
      urljn = 0.d0     
      do j = 1,Nm

c     sA-mA         
      trho2 = rho2(tx,ty,xm(j),ym(j))
      dz = zm(j)-tz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2   
      if (tr2.lt.dms2) goto 712
c     urljn = (dms2/tr2)**6+urljn

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urnew = urnew-terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urnew = -(dielc/tdielc-1.d0)*rtr+urnew
      endif                                                
      endif
      
c      uminew = uminew-rtr

c     sA-mA
      trho2 = rho2(ttx,tty,xm(j),ym(j))
      dz = zm(j)-ttz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2
c     urljo = (dms2/tr2)**6+urljo

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urold = urold-terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urold = -(dielc/tdielc-1.d0)*rtr+urold
      endif                                                
      endif
          
c      umiold = umiold-rtr
      enddo
      urnew = urnew*rnmval
      urold = urold*rnmval
      uminew = uminew*rnmval
c      umiold = umiold*rnmval      

      do j = 1,ip-1
c     sA-sA
      trho2 = rho2(tx,ty,xs(j),ys(j))
      dz = zs(j)-tz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2   
      if (tr2.lt.dss2) goto 712
c     urljn = (dss2/tr2)**6+urljn

      if (tr2.lt.edsq) then           
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urnew = urnew+terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urnew = (dielc/tdielc-1.d0)*rtr+urnew
      endif                                   
      endif      
c      uminew = uminew+rtr

c     sA-sA
      trho2 = rho2(ttx,tty,xs(j),ys(j))
      dz = zs(j)-ttz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2
c     urljo = (dss2/tr2)**6+urljo

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urold = urold+terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urold = (dielc/tdielc-1.d0)*rtr+urold
      endif                                                
      endif
      
c      umiold = umiold+rtr
      enddo

      
      do j = ip+1,Ns
c     sA-sA
      trho2 = rho2(tx,ty,xs(j),ys(j))
      dz = zs(j)-tz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2   
      if (tr2.lt.dss2) goto 712
c     urljn = (dss2/tr2)**6+urljn

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urnew = urnew+terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urnew = (dielc/tdielc-1.d0)*rtr+urnew
      endif                                                
      endif
      
c      uminew = uminew+rtr

c     sA-sA
      trho2 = rho2(ttx,tty,xs(j),ys(j))
      dz = zs(j)-ttz
      dz = dz-h*anint(dz*rh)            
      tr2 = trho2+dz**2
c     urljo = (dss2/tr2)**6+urljo

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      rtr = 1.d0/tr
      urold = urold+terfc(alpha*tr)*rtr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      urold = (dielc/tdielc-1.d0)*rtr+urold
      endif                                                
      endif
      
c      umiold = umiold+rtr
      enddo

      deltaur = bjerrum*(urnew-urold)
c      dumi = bjerrum*(uminew-umiold)      

      deltaukz = 0.d0
      do kntot = 1,knt0
      xrk = vxrk(kntot)   
      yrk = vyrk(kntot)   
      zrk = vzrk(kntot)   
      targ = xrk*tx+yrk*ty+zrk*tz
      ttarg = xrk*ttx+yrk*tty+zrk*ttz
c      rtarg = xrk*tx+yrk*ty+zrk*rtz
c      rttarg = xrk*ttx+yrk*tty+zrk*rttz            
c      dcsum = dcos(ttarg)-dcos(targ)
c      dssum = dsin(ttarg)-dsin(targ)
      dcsum = dcos(ttarg)-dcos(targ)
      dssum = dsin(ttarg)-dsin(targ)
      tvcsum(kntot) = vcsum(kntot)+dcsum
      tvssum(kntot) = vssum(kntot)+dssum
      dc = tvcsum(kntot)**2-vcsum(kntot)**2
      ds = tvssum(kntot)**2-vssum(kntot)**2
      deltaukz = fvec(kntot)*(dc+ds)+deltaukz
      enddo

      deltauk = 0.d0
      do kntot = knt0+1,kntmax
      xrk = vxrk(kntot)   
      yrk = vyrk(kntot)   
      zrk = vzrk(kntot)
      targ = xrk*tx+yrk*ty+zrk*tz
      ttarg = xrk*ttx+yrk*tty+zrk*ttz
c      rtarg = xrk*tx+yrk*ty+zrk*rtz
c      rttarg = xrk*ttx+yrk*tty+zrk*rttz            
c      dcsum = dcos(ttarg)-dcos(targ)+dcos(rtarg)-dcos(rttarg)
c      dssum = dsin(ttarg)-dsin(targ)+dsin(rtarg)-dsin(rttarg)
      dcsum = dcos(ttarg)-dcos(targ)
      dssum = dsin(ttarg)-dsin(targ)
      tvcsum(kntot) = vcsum(kntot)+dcsum
      tvssum(kntot) = vssum(kntot)+dssum
      dc = tvcsum(kntot)**2-vcsum(kntot)**2
      ds = tvssum(kntot)**2-vssum(kntot)**2
      deltauk = fvec(kntot)*(dc+ds)+deltauk
      enddo
      deltauk = 2.d0*deltauk+deltaukz
      deltauk = twopi*deltauk*rvolume
      deltauk = deltauk*bjerrum
      goto 716
 712  idonk = 1
 716  continue
      return
      end      


      subroutine check
      implicit double precision (a-h,o-z)   
      include 'yy.inc'                           
c     k space part:
      zsum = 0.d0
      do kntot = 1,knt0
      trhosq = vcsum(kntot)**2+vssum(kntot)**2
      zsum = zsum+fvec(kntot)*trhosq
      enddo
      sum = 0.d0
      do kntot = knt0+1,kntmax
      trhosq = vcsum(kntot)**2+vssum(kntot)**2
      sum = sum+fvec(kntot)*trhosq
      enddo
      sum = 2.d0*sum+zsum
      sum = twopi*sum*rvolume
      uk = bjerrum*sum
c     real space part:
      ur = 0.d0
      cur = 0.d0
      umi = 0.d0
      urlj = 0.d0
      
      do i = 1,Nm
c     mA:
      tx = xm(i)
      ty = ym(i)
      tz = zm(i)
      if (dabs(tx).gt.ed) write(*,*) '|TX(M)|.GT.ED!!!',i,tx
      if (dabs(ty).gt.ed) write(*,*) '|TY(M)|.GT.ED!!!',i,ty
      if (dabs(tz).gt.halfh) write(*,*) '|TZ(M)|.GT.HALFH!!!',i,tz      
      
      do j = i+1,Nm
c     mA-mA
      tr2 = r2(tx,ty,tz,xm(j),ym(j),zm(j))
      if (tr2.lt.dmm2) write(*,*) 'TR2(MM) LT DMM2',i,j,tr2
c     urlj = (dmm2/tr2)**6+urlj

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      cur = cur+rnval2*terfc(alpha*tr)/tr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      cur = rnval2*(dielc/tdielc-1.d0)/tr+cur
      endif 
      umi = umi+rnval2/tr
      endif
      
      enddo
      
      do j = 1,Ns
c     mA-sA         
      tr2 = r2(tx,ty,tz,xs(j),ys(j),zs(j))         
      if (tr2.lt.dms2) write(*,*) 'TR2(MS) LT DMS2',i,j,tr2
c     urlj = (dms2/tr2)**6+urlj

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      cur = cur-rnmval*terfc(alpha*tr)/tr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      cur = -rnmval*(dielc/tdielc-1.d0)/tr+cur
      endif       
      umi = umi-rnmval/tr
      endif
      
      enddo

      enddo

      do i = 1,Ns
c     sA:
      tx = xs(i)
      ty = ys(i)
      tz = zs(i)
c      rtz= rzs(i)
      if (dabs(tx).gt.ed) write(*,*) '|TX(S)|.GT.ED!!!',i,tx
      if (dabs(ty).gt.ed) write(*,*) '|TY(S)|.GT.ED!!!',i,ty
      if (dabs(tz).gt.halfh) write(*,*) '|TZ(S)|.GT.HALFH!!!',i,tz
      
      do j = i+1,Ns
c     sA-sA:
      tr2 = r2(tx,ty,tz,xs(j),ys(j),zs(j))
      if (tr2.lt.dss2) write(*,*) 'TR2(SS) LT DSS2',i,j,tr2
c     urlj = (dss2/tr2)**6+urlj

      if (tr2.lt.edsq) then            
      tr = dsqrt(tr2)
      cur = cur+terfc(alpha*tr)/tr
      if (tr.lt.rcut) then
      tdielc = cdielc+eratio*(tr-dhspp)
      cur = (dielc/tdielc-1.d0)/tr+cur
      endif       
      umi = umi+1.d0/tr
      endif
      
      enddo

      enddo

      cur = bjerrum*cur
      umi = bjerrum*umi
c      cur = 0.5d0*cur
c      umi = 0.5d0*umi
c      uk = 0.5d0*uk
      ur = cur
      return
      end

      subroutine cinit
      implicit double precision (a-h,o-z)      
      include 'yy.inc'      
      do i = 1,Nm
      tx = (ran2(idum)-0.5d0)*ed2
      ty = (ran2(idum)-0.5d0)*ed2
      tz = (ran2(idum)-0.5d0)*h
c      do j = 1,i-1
c      tr2 = r2(tx,ty,tz,xm(j),ym(j),zm(j))
c      if (tr2.lt.dmm2) goto 636
c      enddo
      xm(i) = tx
      ym(i) = ty
      zm(i) = tz
      enddo
      write(*,*) 'm ions inserted'
      do i = 1,Ns
      tx = (ran2(idum)-0.5d0)*ed2
      ty = (ran2(idum)-0.5d0)*ed2
      tz = (ran2(idum)-0.5d0)*h
c      do j = 1,Nm
c      tr2 = r2(tx,ty,tz,xm(j),ym(j),zm(j))
c      if (tr2.lt.dms2) goto 836
c      enddo
c      do j = 1,i-1
c      tr2 = r2(tx,ty,tz,xs(j),ys(j),zs(j))
c      if (tr2.lt.dss2) goto 836
c      enddo
      xs(i) = tx
      ys(i) = ty
      zs(i) = tz
      enddo
      write(*,*) 's ions inserted'
      return
      end

      function terfc(x)
      implicit double precision (a-h,o-z)   
      include 'yy.inc'                           
      t = 1.d0/(1.d0+0.3275911d0*x)
      terfc = t*(a1+t*(a2+t*(a3+t*(a4+t*a5))))*dexp(-x*x)
      return
      end

      function uc(xx,yy,zz,xxx,yyy,zzz)
      implicit double precision (a-h,o-z)
      include 'yy.inc'
      dx = abs(xx-xxx)
      dx = dx-aint(dx*red)*ed2
      dy = abs(yy-yyy)
      dy = dy-aint(dy*red)*ed2
      dz = abs(zz-zzz)
      dz = dz-aint(dz*red)*ed2
      rmm2 = dx*dx+dy*dy+dz*dz      
      if (rmm2.lt.edsq) then
      r = dsqrt(rmm2)
      uc = terfc(alpha*r)/r
      if (r.lt.rcut) then
      tdielc = cdielc+eratio*(r-dhspp)
      uc = (dielc/tdielc-1.d0)/r+uc
      endif      
      else
      uc = 0.d0
      endif      
      return
      end

      

      function r2(x,y,z,xx,yy,zz)
      implicit double precision (a-h,o-z)   
      include 'yy.inc'                           
      dx = x-xx
      dy = y-yy
      dz = z-zz
      dx = dx-ed2*anint(dx*red2)
      dy = dy-ed2*anint(dy*red2)
      dz = dz-h*anint(dz*rh)      
      r2 = dx*dx+dy*dy+dz*dz
      return
      end

      function dhr2(x,y,z,xx,yy,zz)
      implicit double precision (a-h,o-z)   
      include 'yy.inc'                           
      dx = x-xx
      dy = y-yy
      dz = z-zz
      dx = dx-ed2*anint(dx*red2)
      dy = dy-ed2*anint(dy*red2)
      dz = dz-th*anint(dz*trh)      
      dhr2 = dx*dx+dy*dy+dz*dz
      return
      end      

      function pr2(x,y,z,xx,yy,zz)
      implicit double precision (a-h,o-z)   
      include 'yy.inc'                           
      dx = x-xx
      dy = y-yy
      dz = z-zz
c      dx = dx-ed2*anint(dx*red2)
c      dy = dy-ed2*anint(dy*red2)
c      dz = dz-h*anint(dz*rh)      
      pr2 = dx*dx+dy*dy+dz*dz
      return
      end      


      function rho2(x,y,xx,yy)
      implicit double precision (a-h,o-z)   
      include 'yy.inc'                           
      dx = x-xx
      dy = y-yy
      dx = dx-ed2*anint(dx*red2)
      dy = dy-ed2*anint(dy*red2)
      rho2 = dx*dx+dy*dy
      return
      end

      function rho(x,y,xx,yy)
      implicit double precision (a-h,o-z)   
      include 'yy.inc'                            
      dx = x-xx
      dy = y-yy
      dx = dx-ed2*anint(dx*red2)
      dy = dy-ed2*anint(dy*red2)
      rho = dsqrt(dx*dx+dy*dy)
      return
      end

      function dist(x,y,z,xx,yy,zz)
      implicit double precision (a-h,o-z)   
      include 'yy.inc'                           
      dx = x-xx
      dy = y-yy
      dz = z-zz
      dx = dx-ed2*anint(dx*red2)
      dy = dy-ed2*anint(dy*red2)
      dz = dz-ed2*anint(dz*red2)      
      dist = dsqrt(dx*dx+dy*dy+dz*dz)
      return
      end            
