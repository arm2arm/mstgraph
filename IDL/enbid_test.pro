pro getcicimage, vertx,verty,vertz, wcic, pos, vel,xcic6d

grid=256.0

wcic=vertx*0.0
np=n_elements(wcic)
mx=[min(vertx),max(vertx)]
my=[min(verty),max(verty)]
mz=[min(vertz),max(vertz)]

x=long((vertx-mx[0])/(mx[1]-mx[0])*(grid-1))
y=long((verty-my[0])/(my[1]-my[0])*(grid-1))
z=long((vertz-mz[0])/(mz[1]-mz[0])*(grid-1))

w=wcic+1.0
image=cic(w,x, grid, y, grid, z, grid) 

dV=((mx[1]-mx[0])*(my[1]-my[0])*(mz[1]-mz[0]))/(grid^3.0)
smimage=smooth(image, 6)
for i=0l,np-1 do begin
    wcic[i]=smimage[x[i],y[i],z[i]]
endfor
smimage=0
;return;
;;; Do the 6D cic
gridX=35
gridV=40
mx=[min(pos(0,*)),max(pos(0,*))]
my=[min(pos(1,*)),max(pos(1,*))]
mz=[min(pos(2,*)),max(pos(2,*))]
mvx=[min(vel(0,*)),max(vel(0,*))]
mvy=[min(vel(1,*)),max(vel(1,*))]
mvz=[min(vel(2,*)),max(vel(2,*))]

x=long((pos(0,*)-mx[0])/(mx[1]-mx[0])*(gridX-1))
y=long((pos(1,*)-my[0])/(my[1]-my[0])*(gridX-1))
z=long((pos(2,*)-mz[0])/(mz[1]-mz[0])*(gridX-1))

vx=long((vel(0,*)-mvx[0])/(mvx[1]-mvx[0])*(gridV-1))
vy=long((vel(1,*)-mvy[0])/(mvy[1]-mvy[0])*(gridV-1))
vz=long((vel(2,*)-mvz[0])/(mvz[1]-mvz[0])*(gridV-1))

image6d=fltarr(gridX, gridX, gridX, gridV, gridV, gridV)
for i=0l,np-1 do begin
 ;   print, i, x[i], y[i], z[i]
    image6d[x[i],y[i],z[i], vx[i],vy[i],vz[i]]+=1.0
endfor
xcic6d=wcic*0.0
image6d=smooth(image6d,2, /edge_truncate)
;image6d/=abs((((mx[1]-mx[0])*(my[1]-my[0])*(mz[1]-mz[0])* $
       ;    (mvx[1]-mvx[0])*(mvy[1]-mvy[0])*(mvz[1]-mvz[0]))/gridV^3.0*gridX^3.0))
for i=0l,np-1 do begin
 ;   print, i, x[i], y[i], z[i]
  xcic6d[i]=image6d[x[i],y[i],z[i], vx[i],vy[i],vz[i]]
endfor

end

pro draw_point_image, xvec, yvec, w, xr, yr , ZMINVAL, ZMAXVAL

loadct, 0
plot, xvec, yvec, /nodata, xrange=xr, yrange=yr,  BACKGROUND = 255, COLOR=0
col=alog10(w)
col=(col-alog10(ZMINVAL))/(alog10(ZMAXVAL)-alog10(ZMINVAL))

idz=where(col lt 0, npc)
if(npc gt 0)then COL[idz]=0
idz=where(col gt 1, npc)
if(npc gt 0)then COL[idz]=1


col*=255
idx=sort(col)
np= n_elements(xvec)

loadct, 13
for i=0l,np-1 do begin
ip=idx(i)
oplot, [xvec(ip)], [yvec[ip]], psym=3, color=col(ip)
endfor

end

pro get_new_image, vertx, verty,w,grid, image, xv, yv, nolog=nolog, aver=aver


if(not( IS_DEF(grid))) then grid=256.0

mx=[min(vertx),max(vertx)]
my=[min(verty),max(verty)]

x=(vertx-mx[0])/(mx[1]-mx[0])*0.99*grid
y=(verty-my[0])/(my[1]-my[0])*0.99*grid
image=fltarr(grid, grid)
stop
;image=cic(w,x, grid, y, grid, /average)*1.0e9
if(NOT( IS_DEF(aver)))then begin
    image=cic(w,x, grid, y, grid, /ISOLATED) 
endif else image=cic(double(w),double(x), grid, double(y), grid, /ISOLATED, /AVERAGE)

if(NOT( IS_DEF(nolog)))then begin
    image=image/mean(image,/double)
    image=alog10(image+1)
endif

xv=findgen(grid)/grid *(mx[1]-mx[0])+mx[0]
yv=findgen(grid)/grid *(my[1]-my[0])+my[0]

end


pro plot_two, pos, RR, Ap,  est, vel, ngbp, range=range

vrange=500l
Rplot=5
th=2
color=13
 col='green'
wl=smartlog(est)
xw=findgen(100)/100.0*(max(wl)-min(wl))+min(wl)
h=histogram(wl, nbin=100)
v=max(h,im)

;mi=min(wl)+(xw(im)-min(wl))*0.8
IF KEYWORD_SET(range) THEN BEGIN
temp=smartlog(range)
mi=temp[0]
ma=temp[1]
endif else begin 
  mi=mean(wl)
  ma=max(wl)
endelse
wl=(wl-mi)/(ma-mi)*255
wl=wl>0
iiwl=where(wl gt 255,ng)
if ng gt 0 then wl(iiwl)=255
iwl=sort(wl)


loadct, 0
!P.BACKGROUND=fsc_color('white')
!P.COLOR=fsc_color('black')
plot,pos(0,*), pos(1,*),/xs,/ys,/nodata, psym=3, xrange=[-1,1]*Rplot, yrange=[-1,1]*Rplot
loadct, color
plot_by_points, pos(0,iwl), pos(1,iwl), wl(iwl)
;plots,pos[0,ngbp], pos[1,ngbp], psym=7, thick=th,color=fsc_color(col)

loadct, 0

plot, flatarr(RR), flatarr(Ap),/xs,/ys, psym=3, /nodata, xrange=[0,1]*Rplot
loadct, color
plot_by_points, RR(iwl), Ap(iwl), wl(iwl)
;plots,RR[ngbp], Ap[ngbp], psym=7, thick=th,color=fsc_color(col)

plot, flatarr(vel(0,iwl)), flatarr(vel(1,iwl)),/xs,/ys, psym=3, /nodata, xrange=[-1,1]*vrange, yrange=[-1,1]*vrange
loadct, color
plot_by_points, vel(0,iwl), vel(1,iwl), wl(iwl)
;plots,vel[0,ngbp], vel[1,ngbp], psym=7, thick=th,color=fsc_color(col)
;plots,vel[0,idn1], vel[1,idn1], psym=2, thick=3,color=fsc_color('cyan')
;loadct, color
plot, flatarr(vel(0,iwl)), flatarr(vel(2,iwl)),/xs,/ys, psym=3, /nodata,xrange=[-1,1]*vrange, yrange=[-1,1]*vrange
plot_by_points, vel(0,iwl), vel(2,iwl), wl(iwl)
;return
;plots,vel[0,ngbp], vel[2,ngbp], psym=7, thick=th,color=fsc_color(col)
sc=[0,1,0,1]
sc=[min(RR), max(RR), min(Ap), max(Ap)]
RRs=(RR-sc[0])/(sc[1]-sc[0])
Aps=(Ap-sc[2])/(sc[3]-sc[2])
plot, flatarr(RRs), flatarr(Aps),/xs,/ys, psym=3, /nodata;, xrange=[0,1]*Rplot
loadct, color


openw,1,'c:/arm2arm/DATA/test.u'
nl=n_elements(RRs)
ndim=2l
data=[[RRs],[Aps]]
data[*,0]-=mean(RRs)
data[*,1]-=mean(Aps)
InvCov =  IMSL_INV((Correlate(transpose(data), /COVARIANCE, /DOUBLE))) 
dista=fltarr(nl)
for i=0l,nl-1 do begin
;dista[i]=data[i,*]#InvCov#flatarr(transpose(data[i,*]))*wl[i]
dista[i]=wl[i];*sqrt((data[i,0]^2+data[i,1]^2))
endfor
dista=(dista-min(dista))/(max(dista)-min(dista))*255.0
iso=sort(dista)
plot_by_points, RRs(iso), Aps(iso), dista(iso)
;plots,RRs[idn], Aps[idn], psym=2, thick=3,color=fsc_color('green')
;plots,RRs[idn1], Aps[idn1], psym=2, thick=3,color=fsc_color('cyan')

loadct, 0
plot,pos(0,*), pos(1,*),/xs,/ys,/nodata, psym=3, xrange=[-1,1]*Rplot, yrange=[-1,1]*Rplot
loadct, color
plot_by_points, pos(0,iso), pos(1,iso), dista(iso)

writeu,1,ndim,nl
writeu,1, RRs, Aps
close, 1


;plots,pos[0,idn1], pos[2,idn1], psym=2, thick=3,color=fsc_color('cyan')
;plots,pos[0,600], pos[2,600], psym=2, thick=3,color=fsc_color('blue')

grid=256.0
w=dista
get_image, flatarr(RR), flatarr(Ap),w,grid, image, xv, yv, /nolog, /aver
contour, smooth(image,4), xv, yv, /fill, nlevels=100, MIN_VALUE=mean(image), /xs, /ys

end

pro read_smoothest, file, smest
close, 2
np=0l
openr, 2, file
readu, 2, Np

smest=fltarr(Np)
readu, 2, smest
close, 2

end
pro read_ascii_est,  p,est
base='C:\arm2arm\DATA\MODEL7\MODELS\MODEL7\RUNG2\SNAPS\'
base='/home/arm2arm/Projects/LIA/BAR_FIND/ENBID/'
file=base+'snap_gal_sfr_0450.ascii'
file_est=base+'snap_gal_sfr_0450.ascii_ph4.est'
p=read_ascii(file)
P=P.field1
est=read_ascii(file_est)
est=est.field1
end

pro read_est, file, TYPE,est, ngb, ngbA

close, 1
openr, 1, file
dum1=0l
dum2=0l
make_head, head
readu, 1, dum1, head, dum2
DensNgb=LONG(head.BoxSize)

DensNgbA=LONG(head.time)
;stop
est=fltarr(head.npart(TYPE))
readu, 1, dum1, est,dum2
ngb=lonarr(head.npart(TYPE),DensNgb)
readu, 1, dum1
np=head.npart(TYPE)

ngbt=lonarr(DensNgb)
for in=0l,np-1l do begin
readu, 1, ngbt
ngb(in,*)=ngbt

endfor
return
if(DensNgbA gt 0)then begin

ngba=lonarr(head.npart(TYPE),DensNgbA)
ngbt=lonarr(DensNgbA)
for in=0l,np-1l do begin
;print, in
readu, 1, ngbt
ngbA(in,*)=ngbt

endfor
endif
close, 1

end

pro read_idx, fileidx, grpidx
close, 1
openr, 1, fileidx
ngrp=0l
npID=0l
readu, 1, ngrp
grpidx=INDGEN(ngrp, /ULONG)
grpidx=replicate({Np:0l,pIDS:PTR_NEW()}, ngrp)
for i=0l, ngrp-1 do begin
readu, 1, npID
ids=INDGEN(npID, /ULONG)
readu, 1, ids
grpidx[i].Np=npID
grpidx[i].pIDS=PTR_NEW(ulonarr(npid))
*grpidx[i].pIDS=ids
print, '>>>>>>>', npID
endfor
close, 1

end
pro read_ANNngb, file, annngb
close, 2
np=0l
kl=0l
openr, 2, file
readu, 2,np, kl
annNgb=Replicate({ip:0l,dist:0.0},np)
ip=0l
dist=0.0
readu, 2, annNgb


close, 2


end
pro read_arr, file, arr
 close, 2
 np=0l
 openr, 2, file
 readu, 2, Np

 arr=fltarr(3,Np)
 readu, 2, arr
 close, 2
end

function getWmin, w

return ,10.0^(min(alog10(w))/1.1)
;return ,10.0^(mean(alog10(w)*1.5))
end 

;;;;;;;;;;;;;;;; MAIN ;;;;;;;;;;
device,retain=2,decomposed=0
window,1, xsize=800+400, ysize=800
!P.multi=[0, 6, 4]
Rad=20.0
TYPE=4
mask='snap_gal_sfr'
;mask='snap_gal_0.25_sfr'
;mask='test'
snap='0450'
snap='0950'
base='C:\arm2arm\DATA\MODEL7\MODELS\MODEL7\RUNG2\SNAPS\'
base="/net/direct/dnas01/arm2arm/DATA/LIA/SPHBAR/NFW/MODELS/MODEL7/RUNG2/SNAPS/"
file=base+mask+'_'+snap
file_est=base+mask+'_'+snap+'_ph4.est'
;fileidx='c:\arm2arm\DATA\test.idx'

;read_idx, fileidx, grpidx
snapfile=file
readnew, file, pos, 'POS', parttype=TYPE
readnew, file, vel, 'VEL', parttype=TYPE
readnew, file, m, 'MASS', parttype=TYPE
readnew, file+'_rho_4', rho, 'RHO4'
readnew, file+'_rho_4', hsml, 'H4'
get_com_byPot, snapfile, com

R=Rad
p=pos
p[0,*]-=COM[0]
p[1,*]-=COM[1]
p[2,*]-=COM[2]
d=sqrt(p[0,*]^2.0+p[1,*]^2.0+p[2,*]^2.0)
;ids=sort(d)
idd=where(d lt R, ndd)
print, ndd

hoppath='/net/direct/dnas01/arm2arm/DATA/LIA/SPHBAR/NFW/MODELS/MODEL7/ANALYSIS/HOP/'



file=hoppath+'rho3D64.rho'
read_smoothest, file, smrho


file=hoppath+'rho6D100.rho'
read_smoothest, file, smest

file=hoppath+'rho6D100.rho.sme'
read_smoothest, file, sme

file=hoppath+'rho6D200.rho'
read_smoothest, file, smest200
file=hoppath+'rho6D64.rho'
read_smoothest, file, smest64

file=hoppath+'rho6D64.rho.sme'
read_smoothest, file, smest64sme

file=hoppath+'hsml6D100.hsml'
read_smoothest, file, annhsml

file=hoppath+'mngb6D64.ngb'
read_ANNngb, file, annngb


estfile='/home/arm2arm/Projects/LIA/BAR_FIND/ENBID/snap_gal_sfr_'+snap+'_4_na_myph.est'
read_est, estfile, 4,est, ngb, ngbA




here_test:
cent=COM[0:2];[-9.4345,4.4192,-0.202]

pos(0,*)-=cent[0];mean(pos(0,*), /double)
pos(1,*)-=cent[1];mean(pos(1,*), /double)
pos(2,*)-=cent[2];mean(pos(2,*), /double)

vel(0,*)-=mean(vel(0,idd), /double)
vel(1,*)-=mean(vel(1,idd), /double)
vel(2,*)-=mean(vel(2,idd), /double)
sigma=(vel(0,*)^2.0+vel(1,*)^2.0+vel(2,*)^2.0)/3.0
sigma=sqrt(sigma)

get_velz, pos, vel,  Ar, Ap,Az, RR


nbins=100.0

getcicimage, flatarr(RR),flatarr( Ar), flatarr(Ap),wcic, p, vel,wcic6d
;stop

;; ANN
if(0) then begin
w=smest
Wmin=1e-7
Wmax=0.001
hist=histogram(nbins=nbins,alog10(w), min= alog10(wmin),max=alog10(wmax))
plot,hist,$
  background=fsc_color("white"), color=fsc_color("black"), thick=2


draw_point_image, p[0,*], p[1,*], w, [-1,1]*6, [-1,1]*6 , Wmin,$
  wmax
draw_point_image, vel[0,*], vel[1,*], w, [-1,1]*600, [-1,1]*600 ,Wmin,$
  wmax
draw_point_image, vel[0,*], vel[2,*], w, [-1,1]*600, [-1,1]*600 , Wmin,$
  wmax

draw_point_image, RR[*], Ap[*], w, [0,1]*8, [-1,1]*600 , wmin,$
  wmax

draw_point_image, RR[*], Ar[*], w, [0,1]*8, [-1,1]*600 , wmin,$
  wmax
endif
if(0)then begin
; ANN with smoothing
w=est
Wmin=mean(w)
Wmax=max(w);0.001
hist=histogram(nbins=nbins,alog10(w), min= alog10(wmin),max=alog10(wmax))
plot,hist,$
  background=fsc_color("white"), color=fsc_color("black"), thick=2


draw_point_image, p[0,*], p[1,*], w, [-1,1]*6, [-1,1]*6 , Wmin,$
  wmax
draw_point_image, vel[0,*], vel[1,*], w, [-1,1]*600, [-1,1]*600 ,Wmin,$
  wmax
draw_point_image, vel[0,*], vel[2,*], w, [-1,1]*600, [-1,1]*600 , Wmin,$
  wmax

draw_point_image, RR[*], Ap[*], w, [0,1]*8, [-1,1]*600 , wmin,$
  wmax

draw_point_image, RR[*], Ar[*], w, [0,1]*8, [-1,1]*600 , wmin,$
  wmax
endif
;;;;Primitive CIC
if( 1) then begin
w=wcic

plot, histogram(nbins=nbins,alog10(w), min= alog10(getWmin(w)),max=alog10(max(w))),$
  background=fsc_color("white"), color=fsc_color("black"), thick=2

draw_point_image, p[0,*], p[1,*], w, [-1,1]*Rad, [-1,1]*Rad , getWmin(w),$
  max(w)


draw_point_image, vel[0,*], vel[1,*], w, [-1,1]*600, [-1,1]*600 ,  getWmin(w),$
  max(w)
draw_point_image, vel[0,*], vel[2,*], w, [-1,1]*600, [-1,1]*600 , getWmin(w),$
  max(w)

draw_point_image, RR[*], Ap[*], w, [0,1]*1.5*Rad, [-1,1]*600 , getWmin(w),$
  max(w)

draw_point_image, RR[*], Ar[*], w, [0,1]*1.5*Rad, [-1,1]*600 , getWmin(w),$
  max(w)
endif
if(1)then begin
;; 6D CIC
w=wcic6d

plot, histogram(nbins=nbins,alog10(w), min= alog10(getWmin(w)),max=alog10(max(w))),$
  background=fsc_color("white"), color=fsc_color("black"), thick=2

draw_point_image, p[0,*], p[1,*], w, [-1,1]*Rad, [-1,1]*Rad , getWmin(w),$
  max(w)


draw_point_image, vel[0,*], vel[1,*], w, [-1,1]*600, [-1,1]*600 ,  getWmin(w),$
  max(w)
draw_point_image, vel[0,*], vel[2,*], w, [-1,1]*600, [-1,1]*600 , getWmin(w),$
  max(w)

draw_point_image, RR[*], Ap[*], w, [0,1]*1.5*Rad, [-1,1]*600 , getWmin(w),$
  max(w)

draw_point_image, RR[*], Ar[*], w, [0,1]*1.5*Rad, [-1,1]*600 , getWmin(w),$
  max(w)
endif
if(1) then begin
;; Enbid
w=est
wmin=1e-5
hist=histogram(nbins=nbins,alog10(w), min= alog10(wmin),max=alog10(max(w)))
plot,hist,$
  background=fsc_color("white"), color=fsc_color("black"), thick=2


draw_point_image, p[0,*], p[1,*], w, [-1,1]*Rad, [-1,1]*Rad , Wmin,$
  max(w)
draw_point_image, vel[0,*], vel[1,*], w, [-1,1]*600, [-1,1]*600 ,Wmin,$
  max(w)
draw_point_image, vel[0,*], vel[2,*], w, [-1,1]*600, [-1,1]*600 , Wmin,$
  max(w)

draw_point_image, RR[*], Ap[*], w, [0,1]*Rad*1.5, [-1,1]*600 , wmin,$
  max(w)

draw_point_image, RR[*], Ar[*], w, [0,1]*Rad*1.5, [-1,1]*600 , wmin,$
  max(w)

endif

stop
;; SPH Rho
w=rho
draw_point_image, p[0,*], p[1,*], w, [-1,1]*6, [-1,1]*6 , 10.0^(mean(alog10(w)*1.5)),$
  max(w)



w=smrho
draw_point_image, p[0,*], p[1,*], w, [-1,1]*6, [-1,1]*6 , 10.0^(mean(alog10(w)*1.5)),$
  max(w)

w=sme
draw_point_image, p[0,*], p[1,*], w, [-1,1]*6, [-1,1]*6 , 10.0^(mean(alog10(w)*1.5)),$
  max(w)

w=est
draw_point_image, RR[*], Ap[*], w, [0,1]*8, [-1,1]*600 , 10.0^(mean(alog10(w)*1.5)),$
  max(w)

draw_point_image, RR[*], Ar[*], w, [0,1]*8, [-1,1]*600 , 10.0^(mean(alog10(w)*1.5)),$
  max(w)
;;;; Doing ANNdens test
w=sme
draw_point_image, RR[*], Ap[*], w, [0,1]*8, [-1,1]*600 , 10.0^(mean(alog10(w)*1.5)),$
  max(w)

draw_point_image, RR[*], Ar[*], w, [0,1]*8, [-1,1]*600 , 10.0^(mean(alog10(w)*1.5)),$
  max(w)
w=smest
draw_point_image, RR[*], Ap[*], w, [0,1]*8, [-1,1]*600 , 10.0^(mean(alog10(w)*1.5)),$
  max(w)

draw_point_image, RR[*], Ar[*], w, [0,1]*8, [-1,1]*600 , 10.0^(mean(alog10(w)*1.5)),$
  max(w)


w=smest64
draw_point_image, RR[*], Ap[*], w, [0,1]*8, [-1,1]*600 , 10.0^(mean(alog10(w)*1.5)),$
  max(w)

draw_point_image, RR[*], Ar[*], w, [0,1]*8, [-1,1]*600 , 10.0^(mean(alog10(w)*1.5)),$
  max(w)
w=smest64sme
draw_point_image, RR[*], Ap[*], w, [0,1]*8, [-1,1]*600 , 10.0^(mean(alog10(w)*1.5)),$
  max(w)

draw_point_image, RR[*], Ar[*], w, [0,1]*8, [-1,1]*600 , 10.0^(mean(alog10(w)*1.5)),$
  max(w)

stop








is=sort(RR)
isEst=reverse(sort(est))
Vz=abs(vel(2,*))
smsigma=smooth(sigma(is),10)
;ingb=[155085,158273, 155095]
ip=100l
nn=8
ngbp=0;ngb(is(ip),0:nn)
range=[1e-10,0.0]
;plot_two, pos, RR, Ap,  est,vel, ngbp,range=range
;plot_two, pos, RR, Ap,  smrho/max(smrho),vel, ngbp,range=range
loadct, 39
Grid=256.0
getidinreg, pos, [0,0,0],6, idinreg, sq=1
get_image, flatarr(pos[0,idinreg]), flatarr(pos[1,idinreg]),$
  rho(idinreg),grid, image, xv, yv;, nolog=0, aver=1  

img=smartlog(image)  
contour, img, xv, yv, /fill, nlevels=100,/xs, /ys, MIN_VAL=mean(img) 

get_new_image, flatarr(RR[idinreg]), flatarr(Ap[idinreg]),$
  smest(idinreg)*0.0+1,grid, image, xv, yv;, /nolog,/aver  

imgsm=image;smartlog(image)  
contour, imgsm, xv, yv, /fill, nlevels=100,/xs, /ys, MIN_VAL=mean(imgsm) 
stop
sme=smest/mean(smest)
get_image, flatarr(RR[idinreg]), flatarr(Ap[idinreg]),$
  sme(idinreg),grid, image, xv, yv, nolog=1, aver=1  
stop
img=smartlog(image)  
contour, img, xv, yv, /fill, nlevels=100,/xs, /ys, MIN_VAL=mean(img), MAX_VAL=max(img) 
get_image, flatarr(RR[idinreg]), flatarr(Ar[idinreg]),$
  sme(idinreg),grid, image, xv, yv;, nolog=0, aver=1  
img=smartlog(image)  
contour, img, xv, yv, /fill, nlevels=100,/xs, /ys, MIN_VAL=mean(img), MAX_VAL=max(img) 

stop
;get_image, flatarr(Pos[0,idinreg]), flatarr(Pos[1,idinreg]),$
get_image, flatarr(RR[idinreg]), flatarr(Ap[idinreg]),$
  smest(idinreg),grid, image, xv, yv;,aver=1  ;, nolog=1, aver=1  
imgsm=image;smartlog(image)  
levels=200.0
data=imgsm
step = (Max(data) - Mean(data)) / levels
userLevels = IndGen(levels) * step + Mean(data)
contour, data, xv, yv, /fill, levels=userLevels,C_Colors=Indgen(levels)+3,/xs, /ys;, MIN_VAL=mean(imgsm) 
;plot, histogram(img, nbins=100, min=min(img), max=max(img))
;oplot, histogram(imgsm, nbins=100, min=min(img), max=max(img)), color=fsc_color("blue")
stop
goto, plot_grp
;ngbp=ngba(is(ip),0:nn)
plot_two, pos, RR, Ap,  rho,vel, ngbp,range=range
stop
levels = 200.0
grid=256.0
Rplot=7
loadct, 0
loadct,13
!P.multi=0

;set_teck_eps, 'c:/arm2arm/DATA/all_Enbid.eps', 8, 8
get_image, RR, Ap,est,grid, image, xv, yv;, nolog=nolog, aver=aver
loadct, 0
range=[min(image), max(image)]
step = (range[1] - range[0]) / levels
userLevels = IndGen(levels) * step + range[0]

loadct,4
;reverse_ct
contour, image,xv, yv, /fill, xrange=[0,1]*Rplot, MIN_VALUE=range[0], MAX_VALUE=range[1], /ys, /xs,Levels=userLevels,C_Colors=Indgen(levels)+3

stop
get_image, RR, Ar,est,grid, image, xv, yv;, nolog=nolog, aver=aver
contour, image,xv, yv, /fill, nlevels=200, xrange=[0,1]*Rplot,Levels=userLevels,C_Colors=Indgen(levels)+3
get_image, pos(0,*), pos(1,*),est,grid, image, xv, yv;, nolog=nolog, aver=aver
contour, image,xv, yv, /fill, nlevels=200, xrange=[-1,1]*Rplot, yrange=[-1,1]*Rplot
Rplot=550.
get_image, vel(0,*), vel(1,*),est,grid, image, xv, yv;, nolog=nolog, aver=aver
contour, image,xv, yv, /fill, nlevels=200, xrange=[-1,1]*Rplot, yrange=[-1,1]*Rplot
get_image, vel(0,*), vel(2,*),est,grid, image, xv, yv;, nolog=nolog, aver=aver
contour, image,xv, yv, /fill, nlevels=200, xrange=[-1,1]*Rplot, yrange=[-1,1]*Rplot
get_image, Ap, Ar,est,grid, image, xv, yv;, nolog=nolog, aver=aver
contour, image,xv, yv, /fill, nlevels=200, xrange=[-1,1]*Rplot, yrange=[-1,1]*Rplot
unset_eps

return
hh=histogram(alog10(est), nbin=100)
me=max(hh, i)
lest=alog10(est)
me=findgen(100)/100.*(max(lest)-min(lest))+min(lest)
idm=where(alog10(est)  gt -6)
ngbp=[83,603,836,737]
;plot_two, pos(*,idm), RR(idm), Ap(idm),  est(idm),vel(*, idm), ngbp, range=range
write_png, 'Enbid_all.png', tvrd(true=1)
stop
;plots,RR[ingb], Ap[ingb], psym=7, thick=4 
;plot_two, pos, RR, Ap,  smest/sigma
;plot_two, pos, RR, Ap,  smest/smsigma(is)
;plot_two, pos, RR, Ap,  rho
;plot_two, pos, RR, Ap,  smrho
;plot_two, pos, RR, Ap,  est
is=reverse(is)
plot_grp:
print, grpidx.Np
plot,RR, Ap, psym=7, thick=4,color=fsc_color("black"), xrange=[0,1]*5, yrange=[-1,1]*500, /xs, /ys
loadct,20
col=[ "red", "cyan", "blue", "maroon", "green","black"]
ngrp=n_elements(grpidx)
for i=0l, min([ngrp-1,1]) do begin 
ids=*grpidx[i].pIDS
;;print, ids
;print, i, ids
;oplot, RR(ids), Ap(ids), psym=4, color=fsc_color("red")
oplot, RR(ids), Ap(ids), psym=7, color=fsc_color(col[i mod 5]), thick=2
endfor

;;;;;;;;;;;
;plot,pos[0,*], pos[1,*], psym=7, thick=4,color=fsc_color("black")

loadct,20

ngrp=n_elements(grpidx)
for i=0l, min([ngrp-1,4]) do begin 
ids=*grpidx[i].pIDS
;;print, ids
;print, i, ids
;oplot, RR(ids), Ap(ids), psym=4, color=fsc_color("red")
;oplot, pos[0,ids], pos[1,ids], psym=7, color=fsc_color(col[i mod 5]), thick=2
endfor
;read_arr, 'c:/arm2arm/DATA/ARR.u', ARR
end

