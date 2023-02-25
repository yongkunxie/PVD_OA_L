load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin             

 N_year_beg=1979
 N_year_end=2014
 N_year=N_year_end-N_year_beg+1

 path="~/data"
 
 pi=3.14159
 cp=1004.0
 g=9.8

;---------------------------------------
 filename=path+"/NCEPII/sigma/solar_rad.nc"
 data=addfile(filename,"r")
 sw=data->swhrsig
 dim=dimsizes(sw) 
 print(dim)

 lat=data->Y
 lat=lat/180*pi
 lon=data->X
 lon=lon/180*pi
 sigma=data->Z

 filename=path+"/NCEPII/sigma/surface_pressure.nc"
 data=addfile(filename,"r")
 ps=data->pressfc(0:(dim(0)-1),:,:)

 filename=path+"/NCEPII/sigma/long_rad.nc"
 data=addfile(filename,"r")
 lw=data->lwhrsig

 filename=path+"/NCEPII/sigma/shallow_conv.nc"
 data=addfile(filename,"r")
 sc=data->shahrsig

 filename=path+"/NCEPII/sigma/deep_conv.nc"
 data=addfile(filename,"r")
 dc=data->cnvhrsig

 filename=path+"/NCEPII/sigma/ls_conden.nc"
 data=addfile(filename,"r")
 lc=data->lrghrsig

 filename=path+"/NCEPII/sigma/vert_diff.nc"
 data=addfile(filename,"r")
 vd=data->vdfhrsig

 Q_all=new((/7,dim(0),dim(1),dim(2),dim(3)/),float,-999.0)
 Q_all(0,:,:,:,:)=sw
 Q_all(1,:,:,:,:)=lw
 Q_all(2,:,:,:,:)=sc
 Q_all(3,:,:,:,:)=dc
 Q_all(4,:,:,:,:)=lc
 Q_all(5,:,:,:,:)=vd
 Q_all(6,:,:,:,:)=dim_sum_n(Q_all(0:5,:,:,:,:),0)

 Q_all=Q_all*1.0e6
 do i_term=0,6
  do i_lev=0,dim(1)-1
   Q_all(i_term,:,i_lev,:,:)=Q_all(i_term,:,i_lev,:,:)*ps
  end do
 end do
 Q_all=-Q_all*cp/g
 Q_int=new((/7,dim(0),dim(2),dim(3)/),float,-999.0)
 Q_int=0.0
 do i_lev=0,dim(1)-2
   Q_int=Q_int+(Q_all(:,:,i_lev,:,:)+Q_all(:,:,i_lev+1,:,:))/2.0*(sigma(i_lev+1)-sigma(i_lev))
 end do

 filename_out=path+"/heating/data/NCEPII/Q.vert_int.nc"
 dim1=ispan(1,7,1)
 dim1@units="sw,lw,sc,dc,lc,vd,all"
 dim2=data->T
 copy_VarAtts(data->T,dim2)
 dim4=data->Y
 copy_VarAtts(data->Y,dim4)
 dim5=data->X
 copy_VarAtts(data->X,dim5)
 var=Q_int
 delete(Q_int)
 var@units="W m-2"
;-------------------------------------------------------------------------                    
   system("/bin/rm -f " + filename_out)    ; remove if exists
   fout  = addfile (filename_out, "c")  ; open output file
   setfileoption(fout,"DefineMode",True)
   dimNames = (/"term","time","latitude","longitude"/)  
   dimSizes = (/ dimsizes(dim1) ,dimsizes(dim2) ,dimsizes(dim4) ,dimsizes(dim5)/) 
   dimUnlim = (/ False, False,False,False/)   
   filedimdef(fout,dimNames,dimSizes,dimUnlim)                      
   filevardef(fout, "term"  ,typeof(dim1),dimNames(0))
   filevardef(fout, "time"  ,typeof(dim2),dimNames(1))
   filevardef(fout, "latitude"  ,typeof(dim4),dimNames(2))
   filevardef(fout, "longitude"  ,typeof(dim5),dimNames(3))                
   filevardef(fout, "Q"   ,typeof(var),dimNames)
   filevarattdef(fout,"term"  ,dim1)
   filevarattdef(fout,"time"  ,dim2)
   filevarattdef(fout,"latitude"  ,dim4)
   filevarattdef(fout,"longitude"  ,dim5)  
   filevarattdef(fout,"Q"  ,var)  
   setfileoption(fout,"DefineMode",False)
   fout->term    = (/dim1/) 
   fout->time    = (/dim2/) 
   fout->latitude    = (/dim4/) 
   fout->longitude    = (/dim5/) 
   fout->Q    = (/var/) 
   delete(var)
   delete(dim1)
   delete(dim2)
   delete(dim4)
   delete(dim5)
;------------------------------------------------------------------------- 

end
