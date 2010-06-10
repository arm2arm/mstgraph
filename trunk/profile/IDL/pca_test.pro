pro plot_bh_arrow, epsfile
set_teck_eps, epsfile, 12, 8
!P.thick=6
!P.charthick=4
!P.charsize=2

ff='/dnas01/arm2arm/DATA/LIA/SPHBAR/NFW/MODEL8_BH_EXPLORE/BH/tg1_all.txt'
data1=read_ascii(ff, comment_symbol='#')
data1=data1.field1

ff='/dnas01/arm2arm/DATA/LIA/SPHBAR/NFW/MODEL8_BH_EXPLORE/BH/tg2_all.txt'
data2=read_ascii(ff, comment_symbol='#')
data2=data2.field1
namevec=strarr(16)
cols=fltarr(16)
loadct, 0
plot_io, [0], [0], /xs, /ys, xrange=[0, 6000],  yrange=[0.5, 20], xtitle='Mgas', ytitle='Mbh/Mbhinit'
loadct, 39
init=4
;final=4
for imodel=0, 15 do begin
    cols[imodel]=imodel/15.0*240+10
    ARROW, data1[1,imodel], data1[init,imodel]/data1[3,imodel],$
      data2[1,imodel], data2[init,imodel]/data2[3,imodel],$
      /data, thick=6, /solid, hthick=2, color=cols[imodel]
    namevec[imodel]='ISOLMODEL'+flt2str(imodel)

endfor


old=!p.charthick
!p.charthick=3
!P.charsize=1
tx=2

legend, namevec, /right, color=cols,textcolors=cols
!p.charthick=old
unset_eps

end

pro plot_dump, epsfile, column, yname, yr
set_teck_eps, epsfile, 12, 8
!P.thick=6
!P.charthick=4
!P.charsize=2
namevec=strarr(17)
cols=intarr(17)
loadct, 0
plot , [0], [0], /nodata, xrange=[100,1000], yrange=yr,/xs, /ys, xtitle='Time', Ytitle='N '+yname
loadct, 39
for imodel=0, 15 do begin
LogBase='/dnas01/arm2arm/DATA/LIA/SPHBAR/NFW/MODEL8_BH_EXPLORE/BH/ISOLMODEL'+flt2str(imodel)+'/'
f=logbase+'/dump.log'
print, f
data=read_ascii(f, comment_symbol='#')
data=data.field01

oplot, data[0, *], smooth(data[column, *], 4), color=imodel/15.0*255,LINESTYLE=0

cols[imodel]=imodel/15.0*240+10
namevec[imodel]='ISOLMODEL'+flt2str(imodel)

endfor

old=!p.charthick
!p.charthick=3
!P.charsize=1
tx=2

legend, namevec, /right, color=cols,textcolors=cols
!p.charthick=old
unset_eps
end
pro get_eig, eig , file
eig={eig_str, eigval:fltarr(3), e1:fltarr(3), e2:fltarr(3), e3:fltarr(3)}
data=read_ascii(file)
eig.eigval=data.field01[0:2]
eig.e1=data.field01[3:5]
eig.e2=data.field01[6:8]
eig.e3=data.field01[9:11]
end

function len, a
return, sqrt(a[0]*a[0]+a[1]*a[1])
end

device,retain=2,decomposed=0
;;;;;;;;;;;;;;;;;;;; Windows Models
ImageBase='c:\arm2arm\'
SnapBase='C:\arm2arm\DATA\MODEL7\MODELS\MODEL7\RUNG2\SNAPS\'
LogBase='C:\Documents and Settings\arm2arm\Mes documents\Visual Studio 2008\Projects\mstgraph\profile\'
SnapName='snap_gal_sfr_0450'
;;;;;;;;;;;;;;;;;;;;;
;;Linux Models....
model='04'
ImageBase='./'
SnapName='snap_m8_gal_sfr_bh_450'
SnapBase='/dnas01/arm2arm/DATA/LIA/SPHBAR/NFW/MODEL8_BH_EXPLORE/BH/ISOLMODEL'+model+'/SNAPS/'
LogBase='/dnas01/arm2arm/DATA/LIA/SPHBAR/NFW/MODEL8_BH_EXPLORE/BH/ISOLMODEL'+model+'/'

fAmrad=LogBase+'AmRad.log'
Amrad=read_ascii(fAmrad, comment_symbol='#')

file=SnapBase+snapname
if( FILE_TEST(file))then begin
readnew, file, pos, "POS", parttype=4
if( FILE_TEST(file+'_rho_4')) then $
readnew, file+'_rho_4', rhoin, "RHO4" else rhoin=pos(0,*)*0+1.0
fas=LogBase+'part_ig0.ascii'
data=read_ascii(fas)
data=DOUBLE(data.field1)
COM=[mean(pos[0,*]), mean(pos[1,*]), mean(pos[2,*])]
;COM=[-5.18308, 4.52368, -0.89297]*0
x=flatarr(pos[0,*]-COM[0])
y=flatarr(pos[1,*]-COM[1])
z=flatarr(pos[2,*]-COM[2])

ccol=alog10(rhoin)
col=(ccol-mean(ccol))/(max(ccol)-mean(ccol))*256.0
col=col>0
ig=where(col gt 255, igt)
if(igt gt 0) then col(ig)=255


is=sort(col)
loadct, 0
window, 1, xsize=400, ysize=400
set_teck_eps, ImageBase+'galaxy_marked_m7_0450.eps',  8, 8
;!P.thick=2
;!P.charthick=2
;!P.charsize=2
plot, x, y, psym=3, xrange=[-1,1]*5, yrange=[-1,1]*5,/nodata, /xs, /ys, xtitle='x [kpc]', ytitle='y [kpc]'
loadct,13
for ip=0l, n_elements(col)-1 do begin
    i=is[ip]
    oplot, [x(i)], [y(i)], psym=3, color=col[i]
endfor

;dataAdjust = Transpose([ [x], [y],[z] ])
x=flatarr(data[0,*])
y=flatarr(data[1,*])
z=flatarr(data[2,*])
np=double(n_elements(x))

ip=randomu(10, np*0.1)*np
oplot, data[0,ip], data[1,ip], psym=4,symsize=0.1, color=fsc_color('red')
oplot, [0], [0], psym=1, symsize=100, color=fsc_color('gray'), linestyle=2



dataAdjust = Transpose([ [x], [y],[z] ])
covMatrix = Correlate(dataAdjust, /COVARIANCE, /DOUBLE)
print, "====== COVAR ======="
print, covMatrix
eigenvalues = EIGENQL(covMatrix, EIGENVECTORS=eigenvectors, /DOUBLE)
print, "====== EIGEN VALUES ======="
print, eigenvalues
print, "====== EIGENVECTORS ======="
print, eigenvectors

PhiBar=rad2deg(atan(eigenvectors[1,0],eigenvectors[0,0]))
print, "A/B", eigenvectors[1,0],eigenvectors[0,0]
print, "phi=", PhiBar, -43+180

;featureVector = eigenvectors
;finalData = featureVector ## Transpose(dataAdjust)
;Jvfinal=featureVector ## Transpose(Jv)
fpoint=LogBase+'dum_points.log'
datap=read_ascii(fpoint)
;oplot, datap.field1[0,*], datap.field1[1,*],psym=3

 get_eig, eig , LogBase+'eigendata.log'
 
 scx=max(abs(data[0,ip]))
 scy=max(abs(data[1,ip]))
xx=scale_vector(findgen(11), -1, 1)
!P.thick=8
;oplot, eig.e1[0]*xx*eig.eigval[0], eig.e1[1]*xx*eig.eigval[0], color=fsc_color("magenta")
;oplot, eig.e2[0]*xx*eig.eigval[1], eig.e2[1]*xx*eig.eigval[1], color=fsc_color("magenta")
tvellipse,abs(scx),abs(scy),0,0,PhiBar,thick=10, /MAJOR, /MINOR, /DATA,  color=fsc_color("cyan")



!P.thick=10
mycircle, 0,0,Amrad.field1[1], fsc_color("red")
mycircle, 0,0,Amrad.field1[3], fsc_color("green")
mycircle, 0,0,Amrad.field1[5], fsc_color("blue")
unset_eps
endif

;; Plot the second image
set_teck_eps, ImageBase+'Rm2_m8_0450.eps',  8, 8
fProf=LogBase+'Amf.log'
print, 'Reading: ', fProf
Am2=read_ascii(fProf, comment_symbol='#')

Am=Am2.field1[0:1,*]
loadct, 0
plot, Am[0,*], Am[1, *], xrange=[0,8], yrange=[0,0.8], xtitle="r[kpc]", ytitle="A!dm2!N"
!P.thick=8
oplot, [0,0]+Amrad.field1[1], [0,1], color=fsc_color("red")
oplot, [0,0]+Amrad.field1[3], [0,1], color=fsc_color("green")
oplot, [0,0]+Amrad.field1[5], [0,1], color=fsc_color("blue")

unset_eps

;;;;;;;;;;;;;;;;;;; plot BH/gas part ;;;;;;;;;;;;;;;
plot_bh_arrow, 'BH_init_final.eps'
;;;;;;;;;;;;;;;;;
; plot_dump, 'Gas-all.eps', 1, 'Gas', [0,7000]
; plot_dump, 'Stars-all.eps', 6, 'Stars', [0,20000]


print, "...done..."
end
