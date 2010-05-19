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

ImageBase='c:\arm2arm\'
SnapBase='C:\arm2arm\DATA\MODEL7\MODELS\MODEL7\RUNG2\SNAPS\'
LogBase='C:\Documents and Settings\arm2arm\Mes documents\Visual Studio 2008\Projects\mstgraph\profile\'
file=SnapBase+'snap_gal_sfr_0450'


readnew, file, pos, "POS", parttype=4
readnew, file+'_rho_4', rhoin, "RHO4"
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
oplot, [0], [0], psym=1, symsize=100, color=fsc_color('grey'), linestyle=2



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

fAmrad=LogBase+'AmRad.log'
;!P.thick=5
Amrad=read_ascii(fAmrad, comment_symbol='#')
circle, 0,0,Amrad.field1[1], fsc_color("red")
circle, 0,0,Amrad.field1[3], fsc_color("green")
circle, 0,0,Amrad.field1[5], fsc_color("blue")
unset_eps
;; Plot the second image
set_teck_eps, ImageBase+'Rm2_m7_0450.eps',  8, 8
fProf=LogBase+'Amf.log'
Am2=read_ascii(fProf, comment_symbol='#')

Am=Am2.field1[0:1,*]
loadct, 0
plot, Am[0,*], Am[1, *], xrange=[0,8], yrange=[0,0.8], xtitle="r[kpc]", ytitle="A!dm2!N"
!P.thick=8
oplot, [0,0]+Amrad.field1[1], [0,1], color=fsc_color("red")
oplot, [0,0]+Amrad.field1[3], [0,1], color=fsc_color("green")
oplot, [0,0]+Amrad.field1[5], [0,1], color=fsc_color("blue")

unset_eps


print, "...done..."
end
