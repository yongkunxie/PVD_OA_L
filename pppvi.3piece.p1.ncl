load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin

;----------------------------------prepare data----------------------------------
 N_year=10
 
 a=6371000.0
 omega_c=todouble(7.292*10^(-5))
 pi=3.14159
 R=287.05
 g=9.8
 cp=1004.6
 kappa=R/cp
 cv=cp-R
 
 N_m=10000
 omega_sor=1.6

 N_lat_beg=1
 N_lat_end=53;58
 N_lat=N_lat_end-N_lat_beg+1
 N_lon_beg=0;120
 N_lon_end=239;237
 N_lon=N_lon_end-N_lon_beg+1

 lat_0=45.0/180.0*pi
 f0=2*omega_c*sin(lat_0)
 gp0=100000.0
 p0=100000.0

 acc_cri=0.1/gp0;0.0001
 
 delta_lat=-1.5/180.0*pi
 delta_lon=1.5/180.0*pi
 
 lev_ind=ispan(0,36,1)
 
 path="~/data"
 
 N_season=3

 filename=path+"/pv_diag/ERAI_p/6hr/data/hgt.pf/trend/hgt.trend.season.1979.2010.nc"
 data=addfile(filename,"r")
 gp=data->hgt(N_season,lev_ind,:,:)
 gp=gp*N_year
 dim=dimsizes(gp)
 
 lat=data->latitude
 print((/lat(N_lat_beg),lat(N_lat_end)/))
 lat=lat/180.0*pi
 lon=data->longitude
 print((/lon(N_lon_beg),lon(N_lon_end)/))
 lon=lon/180.0*pi

 z=(data->level(lev_ind))*100

 filename=path+"/pv_diag/ERAI_p/6hr/script/PPPVI.pf/1979.2010/3p/initial/hgt.pppvi.djf.3piece.hb.nc"
 data=addfile(filename,"r")
 gp_hb=data->hgt(:,lev_ind,:,:)
 gp_hb=gp_hb*N_year
 
 filename=path+"/pv_diag/ERAI_p/6hr/script/PPPVI.pf/1979.2010/3p/initial/hgt.pppvi.djf.3piece.lb.nc"
 data=addfile(filename,"r")
 gp_lb=data->hgt(:,lev_ind,:,:)
 gp_lb=gp_lb*N_year
 
 gp_g_z=new(dim,double,-999.0)
 gp_g_z(1:(dim(0)-2),:,:)=(gp(2:(dim(0)-1),:,:)-gp(0:(dim(0)-3),:,:))/conform_dims((/dim(0)-2,dim(1),dim(2)/),z(2:(dim(0)-1))-z(0:(dim(0)-3)),0)
  
 gp_2g_z=new(dim,double,-999.0)
 gp_2g_z(1:(dim(0)-2),:,:)=p0*(gp(2:(dim(0)-1),:,:)+gp(0:(dim(0)-3),:,:)-2.0*gp(1:(dim(0)-2),:,:))/conform_dims((/dim(0)-2,dim(1),dim(2)/),(z(2:(dim(0)-1))-z(1:(dim(0)-2)))*(z(1:(dim(0)-2))-z(0:(dim(0)-3))),0)
    
 filename=path+"/pv_diag/ERAI_p/6hr/data/ppv.pf/trend/ppv.trend.season.1979.2010.nc"
 data=addfile(filename,"r")
 pv=data->pv(N_season,:,:,:)
 pv=pv*N_year
  
 filename=path+"/pv_diag/ERAI_p/6hr/data/ppv.pf/ppv.ref.S.nc"
 data=addfile(filename,"r")
 rec_s_clim=data->pv(0,3,:)
 rec_s_g_z_clim=data->pv(1,3,:)

;------------------------------------Mask-----------------------------------------
 filename=path+"/pv_diag/ERAI_p/6hr/data/hgt/trend/hgt.trend.season.nc"
 data=addfile(filename,"r")
 mask_ob=data->hgt(N_season,lev_ind,:,:)
 mask_ind=new((/dim(1),dim(2)/),integer,999)
 do i_lat=0,dim(1)-1
  do i_lon=0,dim(2)-1
   do i_lev=0,dim(0)-1
    if (.not. ismissing(mask_ob(i_lev,i_lat,i_lon))) then
     mask_ind(i_lat,i_lon)=i_lev
    end if
   end do
  end do
 end do
;------------------------------------Mask-----------------------------------------
       
;----------------------------nondimensionalization--------------------------------
 gp=gp/gp0
 gp_lb=gp_lb/gp0
 gp_hb=gp_hb/gp0
 pv=pv/(f0*10000.0)
 z=z/p0
 c1=gp0/(a*a*f0*f0)
 c2=rec_s_g_z_clim
 c3=rec_s_clim/p0
;----------------------------nondimensionalization--------------------------------

;----------------------------------prepare data----------------------------------

;-----------------------------------piecewise------------------------------------
 N_piece=1
 pv_temp=pv
 delete(pv)
 pv=new((/N_piece,dim(0),dim(1),dim(2)/),double,-999.0)
 do i_piece=0,N_piece-1
  pv(i_piece,:,:,:)=pv_temp
 end do
 pv(0,0:(dim(0)-3),:,:)=0.0
 delete(pv_temp)
 
 gp_temp=gp
 delete(gp)
 gp=new((/N_piece,dim(0),dim(1),dim(2)/),double,-999.0)
 gp=0.0
 do i_lon=0,dim(2)-1
  do i_piece=0,N_piece-1
   gp(i_piece,:,N_lat_beg-1,i_lon)=gp_lb(0,:,N_lat_beg-1,i_lon)
   gp(i_piece,:,N_lat_end+1,i_lon)=gp_lb(0,:,N_lat_end+1,i_lon)
  end do
 end do
 delete(gp_temp)
 
 rho_temp=-gp_g_z
 rho=new((/N_piece,dim(0),dim(1),dim(2)/),double,-999.0)
 do i_piece=0,N_piece-1
  rho(i_piece,:,:,:)=rho_temp
 end do
 rho(0,0:15,:,:)=0.0
 delete(rho_temp)
 
 rho_g_z_temp=-gp_2g_z
 rho_g_z=new((/N_piece,dim(0),dim(1),dim(2)/),double,-999.0)
 do i_piece=0,N_piece-1
  rho_g_z(i_piece,:,:,:)=rho_g_z_temp
 end do
 rho_g_z(0,0:15,:,:)=0.0
 delete(rho_g_z_temp)
;-----------------------------------piecewise------------------------------------

;-----------------------------------piecewise------------------------------------

 z_k=new(dimsizes(z),double,-999.0)
 z_k(1:(dim(0)-2))=1.0/(z(2:(dim(0)-1))-z(0:(dim(0)-3)))

 z_k_k=new(dimsizes(z),double,-999.0)
 z_k_k(1:(dim(0)-2))=1.0/(z(2:(dim(0)-1))-z(1:(dim(0)-2)))/(z(1:(dim(0)-2))-z(0:(dim(0)-3)))

 c4=c1/(((cos(lat))^2)*(delta_lon^2))
 c5=c1*(1.0/(delta_lat^2)-tan(lat)/(2*delta_lat))
 c6=c1*(1.0/(delta_lat^2)+tan(lat)/(2*delta_lat))
 c7=c2*z_k+c3*z_k_k
 c8=c3*z_k_k-c2*z_k
 c11=2.0*c4+2.0*c1/(delta_lat^2)
 c9=conform_dims(dim,c11,1)+conform_dims(dim,2.0*c3*z_k_k,0)
;-----------------------------------iteration------------------------------------
 N_lev_beg=0
 do i_piece=0,N_piece-1
  
  gp_diff=(/gp(i_piece,:,:,:)+10000.0/)
  
  i_m=1
  do while ((i_m .le. N_m) .and. (max(abs(gp_diff)) .gt. acc_cri))

   gp_diff=0.0

   gp_temp1=gp(i_piece,:,:,:)

  ;*****************************************begin longitude**********************************
   do i_lat=N_lat_beg,N_lat_end
    do i_lon=0,0
    
     N_lev_end=mask_ind(i_lat,i_lon)
     ;----------------------------------upper boundary-------------------------------------
     i_lev=N_lev_beg
     gp(i_piece,i_lev,i_lat,i_lon)=gp_hb(0,i_lev,i_lat,i_lon)
     ;----------------------------------upper boundary-------------------------------------
      
     ;--------------------------------bottom boundary--------------------------------------
     i_lev=N_lev_end
     gp(i_piece,i_lev,i_lat,i_lon)=gp_hb(0,i_lev,i_lat,i_lon)
     ;--------------------------------bottom boundary--------------------------------------

     ;-------------------------------interior boundary-------------------------------------
     do i_lev=N_lev_beg+1,N_lev_end-1

      gp(i_piece,i_lev,i_lat,i_lon)=(1-omega_sor)*gp(i_piece,i_lev,i_lat,i_lon)+omega_sor*(-pv(i_piece,i_lev,i_lat,i_lon)+c4(i_lat)*(gp(i_piece,i_lev,i_lat,i_lon+1)+gp(i_piece,i_lev,i_lat,dim(2)-1))+c5(i_lat)*gp(i_piece,i_lev,i_lat+1,i_lon)+c6(i_lat)*gp(i_piece,i_lev,i_lat-1,i_lon)+c7(i_lev)*gp(i_piece,i_lev+1,i_lat,i_lon)+c8(i_lev)*gp(i_piece,i_lev-1,i_lat,i_lon))/c9(i_lev,i_lat,i_lon)

     end do
     ;-------------------------------interior boundary-------------------------------------
     
     gp_diff(N_lev_beg:N_lev_end,i_lat,i_lon)=(/gp(i_piece,N_lev_beg:N_lev_end,i_lat,i_lon)-gp_temp1(N_lev_beg:N_lev_end,i_lat,i_lon)/)
    
    end do
   end do
  ;*****************************************begin longitude**********************************

  ;******************************************end longitude***********************************
   do i_lat=N_lat_beg,N_lat_end
    do i_lon=dim(2)-1,dim(2)-1
     
     N_lev_end=mask_ind(i_lat,i_lon)
     ;----------------------------------upper boundary-------------------------------------
     i_lev=N_lev_beg
     gp(i_piece,i_lev,i_lat,i_lon)=gp_hb(0,i_lev,i_lat,i_lon)
     ;----------------------------------upper boundary-------------------------------------
      
     ;--------------------------------bottom boundary--------------------------------------
     i_lev=N_lev_end
     gp(i_piece,i_lev,i_lat,i_lon)=gp_hb(0,i_lev,i_lat,i_lon)
     ;--------------------------------bottom boundary--------------------------------------
    
     ;-------------------------------interior boundary-------------------------------------
     do i_lev=N_lev_beg+1,N_lev_end-1

      gp(i_piece,i_lev,i_lat,i_lon)=(1-omega_sor)*gp(i_piece,i_lev,i_lat,i_lon)+omega_sor*(-pv(i_piece,i_lev,i_lat,i_lon)+c4(i_lat)*(gp(i_piece,i_lev,i_lat,0)+gp(i_piece,i_lev,i_lat,i_lon-1))+c5(i_lat)*gp(i_piece,i_lev,i_lat+1,i_lon)+c6(i_lat)*gp(i_piece,i_lev,i_lat-1,i_lon)+c7(i_lev)*gp(i_piece,i_lev+1,i_lat,i_lon)+c8(i_lev)*gp(i_piece,i_lev-1,i_lat,i_lon))/c9(i_lev,i_lat,i_lon)

     end do
     ;-------------------------------interior boundary-------------------------------------
     
     gp_diff(N_lev_beg:N_lev_end,i_lat,i_lon)=(/gp(i_piece,N_lev_beg:N_lev_end,i_lat,i_lon)-gp_temp1(N_lev_beg:N_lev_end,i_lat,i_lon)/)
    
    end do
   end do
  ;******************************************end longitude***********************************

  ;***************************************interior longitude*********************************
   do i_lat=N_lat_beg,N_lat_end
    do i_lon=1,dim(2)-2
    
     N_lev_end=mask_ind(i_lat,i_lon)
     ;----------------------------------upper boundary-------------------------------------
     i_lev=N_lev_beg
     gp(i_piece,i_lev,i_lat,i_lon)=gp_hb(0,i_lev,i_lat,i_lon)
     ;----------------------------------upper boundary-------------------------------------
      
     ;--------------------------------bottom boundary--------------------------------------
     i_lev=N_lev_end
     gp(i_piece,i_lev,i_lat,i_lon)=gp_hb(0,i_lev,i_lat,i_lon)
     ;--------------------------------bottom boundary--------------------------------------
    
     ;-------------------------------interior boundary-------------------------------------
     do i_lev=N_lev_beg+1,N_lev_end-1

      gp(i_piece,i_lev,i_lat,i_lon)=(1-omega_sor)*gp(i_piece,i_lev,i_lat,i_lon)+omega_sor*(-pv(i_piece,i_lev,i_lat,i_lon)+c4(i_lat)*(gp(i_piece,i_lev,i_lat,i_lon+1)+gp(i_piece,i_lev,i_lat,i_lon-1))+c5(i_lat)*gp(i_piece,i_lev,i_lat+1,i_lon)+c6(i_lat)*gp(i_piece,i_lev,i_lat-1,i_lon)+c7(i_lev)*gp(i_piece,i_lev+1,i_lat,i_lon)+c8(i_lev)*gp(i_piece,i_lev-1,i_lat,i_lon))/c9(i_lev,i_lat,i_lon)

     end do
     ;-------------------------------interior boundary-------------------------------------
     
     gp_diff(N_lev_beg:N_lev_end,i_lat,i_lon)=(/gp(i_piece,N_lev_beg:N_lev_end,i_lat,i_lon)-gp_temp1(N_lev_beg:N_lev_end,i_lat,i_lon)/)
    
    end do
   end do
  ;***************************************interior longitude*********************************

   print(max(abs(gp_diff)))
   i_m=i_m+1
  end do
 end do
;-----------------------------------iteration------------------------------------

;-------------------------------restore dimension---------------------------------
 gp=gp*gp0/N_year
;-------------------------------restore dimension---------------------------------

 filename=path+"/pv_diag/ERAI_p/6hr/data/hgt.pf/trend/hgt.clim.season.nc"
 data=addfile(filename,"r")
 
 filename_out=path+"/pv_diag/ERAI_p/6hr/data/PPPVI.pf/1979.2010/3p/hgt.pppvi.djf.3piece.p1.nc"
 dim1=ispan(1,N_piece,1)
 dim1@units="piece numbers"
 dim3=tofloat(data->level(lev_ind))
 copy_VarAtts(data->level,dim3)
 dim4=tofloat(data->latitude)
 copy_VarAtts(data->latitude,dim4)
 dim5=tofloat(data->longitude)
 copy_VarAtts(data->longitude,dim5)
 var=gp
 delete(gp)
 var@units="10-6 PVU"
;-------------------------------------------------------------------------                    
   system("/bin/rm -f " + filename_out)    ; remove if exists
   fout  = addfile (filename_out, "c")  ; open output file
   setfileoption(fout,"DefineMode",True)
   dimNames = (/"piece","level","latitude", "longitude"/)  
   dimSizes = (/ dimsizes(dim1),dimsizes(dim3),dimsizes(dim4),dimsizes(dim5)/) 
   dimUnlim = (/ False,False,False,False/)   
   filedimdef(fout,dimNames,dimSizes,dimUnlim)                      
   filevardef(fout, "piece"  ,typeof(dim1),dimNames(0))                
   filevardef(fout, "level"  ,typeof(dim3),dimNames(1))                
   filevardef(fout, "latitude" ,typeof(dim4),dimNames(2))
   filevardef(fout, "longitude" ,typeof(dim5),dimNames(3))
   filevardef(fout, "hgt"   ,typeof(var),dimNames)
   filevarattdef(fout,"piece"  ,dim1)
   filevarattdef(fout,"level"  ,dim3)
   filevarattdef(fout,"latitude"  ,dim4)  
   filevarattdef(fout,"longitude"  ,dim5)  
   filevarattdef(fout,"hgt"  ,var)  
   setfileoption(fout,"DefineMode",False)
   fout->piece    = (/dim1/) 
   fout->level    = (/dim3/) 
   fout->latitude    = (/dim4/) 
   fout->longitude    = (/dim5/) 
   fout->hgt    = (/var/) 
   delete(var)
   delete(dim1)
   delete(dim3)
   delete(dim4)
   delete(dim5)
;-------------------------------------------------------------------------  

end

