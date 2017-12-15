program one
    implicit none
	real,parameter::G=101.94,p0_0=101325,t0_0=288.15,n=7500,pi_c=4.1,z=8         !初始数据
	real,parameter::k=1.4,k1=0.0404,pi=3.1415926,r=287.06
	!导流叶片
	real::sigma_in=0.985
	real::p1_0,t1_0
	!第一级平均直径上气流及几何参数计算
	real::c1_a1=160,d1_1=0.4                  !自取
	real::p1_1,t1_1,rho1_1,c1_1,q1_1,lambda1_1,d1_t1,d1_h1,a1_1,L1_1  !求解
	!末级平均直径上气流及几何参数计算
	real::eta_c=0.872,c_ac=115,
	real::pc_0,tc_0,pc_1,tc_1,rhoc_1,cc_1,qc_1,lambdac_1,dc_t1,dc_h1,ac_1,Lc_1
	!选择级平均加工量
	 real::h_ad,eta_ad,h_adm,eta_adm
	 real::a_0=0.02332,a_2
	 real::a_1
	 !多级计算过程
	 real,dimension(8)::h_adi,c_a2_a,h_i,eta_a
	 real,dimension(9)::p_01_a,t_01_a,rhoc_11_a,p_11_a,t_11_a,p_02_a,t_02_a,t_12_a,c_a1_a,alpha_1_a,&
	                     &c_1_a,lambda_1_a,q_1_a,s_1,d_m1,d_t1,d_h1,L1
	real,dimension(8)::u_1_a,u_2_a,c_1u_a,c_2u_a,w_1u_a,w_2u_a,w_1_a,w_2_a,c_2_a,alpha_2_a,beta_1_a,&
	                   &beta_2_a,M_w_a,d_m2,delta_t_a,lambda_2_a,rhoc_12_a,s_2,d_t2,d_h2,d_m20,p_12_a,L2
	real,dimension(8)::p_03_a,t_03_a,M_c_a,omega_a,lambda_1w_a,t_01_w_a,p_01_w_a,t_02_w_a,p_02_w_a,p_02_w_a_ad,&
	                   &lambda_2w_a,omega_sunshi_ja,omega_sunshi_da
	 real::c_p
	 integer::i,j
	!确定各列叶栅气动计算尺寸
	real::K_G,delta_h,delta_t
	real,dimension(8)::r_he,r_te,r_1_h,r_2_h,r_1_t,r_2_t,r_1_m,r_2_m
	!确定径向计算站数目及其位置
	real,dimension(7,8)::r_r1,u_1,c_1u,c_1a,c_1,w_1u,w_1,alpha_1,beta_1,lambda_1,t_11,M_w,delta_A_1,delta_G_1,delta_GA_1,&
	                     &p_01,p_11,rhoc_11
	real::z_r
	real,dimension(8)::A,B,eps_1,G_c_1,t_0average_1,p_0average_1,t_0zong_1,p_0zong_1
	real,dimension(7)::sigma_1
	!!动叶出口截面
	real,dimension(8)::E,h,C,D
	real,dimension(7,8)::r_r2,u_2,c_2u,c_2a,c_2,w_2u,w_2
	!计算气体状态（动叶）
	real,dimension(8)::delta_t0,t_02,G_c_2,t_0average_2,p_0average_2,eta_average,t_0zong_2,p_0zong_2,eta_zong
	real,dimension(7,8)::lambda_2,t_12,M_c,t_01_w,t_02_w,lambda_1w,lambda_2w,p_01_w,p_02_w,p_02_w_ad,p_12,p_02,rhoc_12,&
	                        &delta_GA_2,delta_A_2,delta_G_2,omega,eta_adi,alpha_2,beta_2
    real,dimension(7)::eps_2
  
	!导流叶片
	open(unit=10,file='result.txt')
	p1_0=p0_0*sigma_in
    t1_0=t0_0
    write(10,*)'导流叶片:','p1_0=',p1_0,'t1_0=',t1_0
	!第一级平均直径上气流及几何参数计算
	!call jmcs(c1_a1,d1_1,alpha1_1,t1_0,p1_0,p1_1,t1_1,rho1_1,c1_1,q1_1,lambda1_1,d1_t1,d1_h1,a1_1,L1_1)
	c1_1=c1_a1/sin(alpha1_1*pi/180)
	lambda1_1=c1_1/sqrt(2*r*k*t1_0/(k+1))
	t1_1=(1-(k-1)*lambda1_1**2/(k+1))*t1_0
	p1_1=(1-(k-1)*lambda1_1**2/(k+1))**(k/(k-1))*p1_0
	q1_1=lambda1_1*((k+1)/2-(k-1)*(lambda1_1**2)/2)**(1/(k-1))
	a1_1=G*sqrt(t1_0)/(k1*q1_1*p1_0*sin(alpha1_1*pi/180))
	d1_t1=sqrt(4*a1_1/(pi*(1-d1_1**2)))
	d1_h1=d1_t1*d1_1
	L1_1=(d1_t1-d1_h1)/2
	write(10,*)'第一级平均直径上气流及几何参数计算','c1_1=',c1_1,'lambda1_1=',lambda1_1,'t1_1=',t1_1,'p1_1=',p1_1,'q1_1=',q1_1,&
	&'a1_1=',a1_1,'d1_t1=',d1_t1,'d1_h1=',d1_h1,'L1_1=',L1_1
	!末级平均直径上气流及几何参数计算
	pc_0=p1_0*pi_c
	tc_0=t1_0*(1+(pi_c**((k-1)/k)-1)/eta_c)
	cc_1=c_ac
	lambdac_1=cc_1/sqrt(2*r*k*tc_0/(k+1))
	tc_1=(1-(k-1)*lambdac_1**2/(k+1))*tc_0
	pc_1=(1-(k-1)*lambdac_1**2/(k+1))**(k/(k-1))*pc_0
	qc_1=lambdac_1*((k+1)/2-(k-1)*(lambdac_1**2)/2)**(1/(k-1))
	ac_1=G*sqrt(tc_0)/(k1*qc_1*pc_0)
	dc_t1=d1_t1
	dc_h1=sqrt(dc_t1**2-4*ac_1/pi)
	Lc_1=(dc_t1-dc_h1)/2
	write(10,*)'末级平均直径上气流及几何参数计算','cc_1=',cc_1,'lambdac_1=',lambdac_1,'tc_1=',tc_1,'pc_1=',pc_1,'qc_1=',qc_1&
	&,'ac_1=',ac_1,'dc_t1=',dc_t1,'dc_h1=',dc_h1,'Lc_1=',Lc_1,'tc_0=',tc_0,'pc_0=',pc_0
	!选择级平均加工量
	do while(.true.)                                       !迭代求重热系数
	h_ad=k*r*t1_1*((pc_1/p1_1)**((k-1)/k)-1)/(k-1)
	eta_ad=h_ad/(k*r*(tc_1-t1_1)/(k-1))
	eta_adm=eta_ad*(1+a_0)
	a_1=eta_adm*((pc_1/p1_1)**((k-1)/(k*eta_adm))-1)/((pc_1/p1_1)**((k-1)/k)-1)-1
	a_2=a_1*(z-1)/z
	if((abs(a_2-a_0))/a_0<0.0001)exit
	a_0=a_2
	end do
	h_adm=(1+a_2)*h_ad/z
	write(10,*)'选择级平均加工量','h_ad=',h_ad,'eta_ad=',eta_ad,'h_adm=',h_adm,'eta_adm=',eta_adm,'a_2=',a_2
	!多级计算过程
	h_adi(1)=0.524*h_adm                                                        !分配能量
	h_adi(2)=0.8428*h_adm
	h_adi(3)=1.1502*h_adm
	h_adi(4)=1.1729*h_adm
	h_adi(5)=1.1824*h_adm
	h_adi(6)=1.1721*h_adm
	h_adi(7)=0.9883*h_adm
	h_adi(8)=0.9673*h_adm
	c_a1_a=(/160,160,160,160,155,150,145,140,135/)
	c_a2_a=(/158,156,154,152,152,148,143,138/)
	alpha_1_a=(/60.42,63.62,65.26,66.72,68.92,69.87,71.59,73.95,90/)
	eta_a=(/0.86,0.88,0.91,0.908,0.906,0.904,0.892,0.87/)
	c_p=k*r/(k-1)
	t_01_a(1)=t1_0                 !将之前求解的第一级参数赋值
	p_01_a(1)=p1_0
	do i=1,9
	if(i<9)then
	h_i(i)=h_adi(i)/eta_a(i)
	t_01_a(i+1)=h_i(i)/c_p+t_01_a(i)
	p_01_a(i+1)=p_01_a(i)*(h_adi(i)/(c_p*t_01_a(i))+1)**(k/(k-1))
	end if
	c_1_a(i)=c_a1_a(i)/sin(alpha_1_a(i)*pi/180)
	lambda_1_a(i)=c_1_a(i)/sqrt(2*r*k*t_01_a(i)/(k+1))
	t_11_a(i)=(1-(k-1)*lambda_1_a(i)**2/(k+1))*t_01_a(i)
	p_11_a(i)=(1-(k-1)*lambda_1_a(i)**2/(k+1))**(k/(k-1))*p_01_a(i)
	rhoc_11_a(i)=p_11_a(i)/(t_11_a(i)*r)
	q_1_a(i)=lambda_1_a(i)*((k+1)/2-(k-1)*(lambda_1_a(i)**2)/2)**(1/(k-1))
	s_1(i)=G*sqrt(t_01_a(i))/(k1*q_1_a(i)*p_01_a(i)*sin(alpha_1_a(i)*pi/180))
	d_t1(i)=d1_t1                                                       ! d_m1,d_t1,d_h1,L1
	d_h1(i)=sqrt(d_t1(i)**2-4*s_1(i)/pi)
	d_m1(i)=(d_h1(i)+d_t1(i))/2
	L1(i)=(d_t1(i)-d_h1(i))/2
	write(10,*)'进口截面的状态和几何参数','i=',i,'c_1_a=',c_1_a(i),'alpha_1_a=',alpha_1_a(i),'lambda_1_a=',lambda_1_a(i),'t_11_a=',t_11_a(i),&
	&'p_11_a=',p_11_a(i),'q_1_a=',q_1_a(i),'s_1=',s_1(i),'d_t1=',d_t1(i),'d_h1=',d_h1(i),'d_m1=',d_m1(i),'L1=',L1(i)
	end do
	do i=1,8
	c_1u_a(i)=sqrt((c_1_a(i))**2-(c_a1_a(i))**2)
	u_1_a(i)=n*pi*d_m1(i)/60
	w_1u_a(i)=u_1_a(i)-c_1u_a(i)
	w_1_a(i)=sqrt((w_1u_a(i))**2+(c_a1_a(i))**2)
	M_w_a(i)=w_1_a(i)/sqrt(k*r*t_11_a(i))
	if(M_w_a(i)>0.9)then
	write(10,*)"this is wrong"
	end if
	beta_1_a(i)=asin(c_a1_a(i)/w_1_a(i))*180/pi
	d_m20(i)=(d_m1(i)+d_m1(i+1))/2
	do while(.true.)
	u_2_a(i)=n*pi*d_m20(i)/60
	c_2u_a(i)=(h_i(i)+c_1u_a(i)*u_1_a(i))/u_2_a(i)
	c_2_a(i)=sqrt(c_2u_a(i)**2+c_a2_a(i)**2)
	w_2u_a(i)=u_2_a(i)-c_2u_a(i)
	w_2_a(i)=sqrt((w_2u_a(i))**2+(c_a2_a(i))**2)
	beta_2_a(i)=asin(c_a2_a(i)/w_2_a(i))*180/pi
	alpha_2_a(i)=asin(c_a2_a(i)/c_2_a(i))*180/pi
	delta_t_a(i)=(h_i(i)-0.5*(c_2_a(i)**2-c_1_a(i)**2))*(k-1)/(k*r)
	t_12_a(i)=t_11_a(i)+delta_t_a(i)
    p_12_a(i)=p_11_a(i)*(t_12_a(i)/t_11_a(i))**(k*eta_a(i)/(k-1))
	rhoc_12_a(i)=p_12_a(i)/(t_12_a(i)*r)
	s_2(i)=G/(c_a2_a(i)*rhoc_12_a(i))
	d_t2(i)=d1_t1                                                       ! d_m2,d_t2,d_h2,L1
	d_h2(i)=sqrt(d_t2(i)**2-4*s_2(i)/pi)
	d_m2(i)=(d_h2(i)+d_t2(i))/2
	L2(i)=(d_t2(i)-d_h2(i))/2

	if(abs((d_m2(i)-d_m20(i))/d_m20(i))<0.001)then                      !可能进入无限循环
	exit
	end if
	d_m20(i)=d_m2(i)
	end do
	!write(10,*)'动叶出口计算','i=',i,'c_1u_a=',c_1u_a(i),'u_1_a',u_1_a(i),'w_1u_a=',w_1u_a(i),'w_1_a=',w_1_a(i),&
	!&'M_w_a=',M_w_a(i),'beta_1_a=',beta_1_a(i),'u_2_a',u_2_a(i),'c_2u_a=',c_2u_a(i),'c_2_a=',c_2_a(i),'w_2u_a=',w_2u_a(i),&
	!&'w_2_a=',w_2_a(i),'beta_2_a=',beta_2_a(i),'alpha_2_a=',alpha_2_a(i),'delta_t_a=',delta_t_a(i),'t_12_a=',t_12_a(i),&
	!&'p_12_a=',p_12_a(i),'rhoc_12_a=',rhoc_12_a(i),'s_2=',s_2(i),'d_t2=',d_t2(i),'d_h2=',d_h2(i),'d_m2=',d_m2(i)
    
	!write(*,*)d_m20(i),d_m2(i),d_m1(i)
	t_02_a(i)=t_01_a(i)+h_i(i)/c_p
	lambda_2_a(i)=c_2_a(i)/sqrt(2*k*r*t_02_a(i)/(k+1))
	p_02_a(i)=p_12_a(i)/(1-(k-1)*lambda_2_a(i)**2/(k+1))**(k/(k-1))
	p_03_a(i)=p_01_a(i)*(1+(k-1)*h_adi(i)/(k*r*p_01_a(i)))**(k/(k-1))
	t_03_a(i)=t_02_a(i)
	M_c_a(i)=c_2_a(i)/sqrt(k*r*t_12_a(i))
	if(M_c_a(i)>0.9)then
	write(10,*)'c is wrong'
	end if
	omega_a(i)=1-(c_2_a(i)**2-c_1_a(i)**2)/(2*h_i(i))
	omega_sunshi_ja(i)=(p_02_a(i)-p_03_a(i))/(p_02_a(i)-p_12_a(i))
	lambda_1w_a(i)=w_1_a(i)/sqrt(2*k*r*t_01_a(i)/(k+1))
	t_01_w_a(i)=t_11_a(i)/(1-(k-1)*lambda_1w_a(i)**2/(k+1))
	p_01_w_a(i)=p_01_a(i)/(1-(k-1)*lambda_1w_a(i)**2/(k+1))**(k/(k-1))
	t_02_w_a(i)=t_01_w_a(i)*(1+u_2_a(i)**2*(1-(d_m1(i)/d_m2(i))**2)*(k-1)/(2*k*r*t_01_w_a(i)))
	p_02_w_a_ad(i)=p_01_w_a(i)*(t_02_w_a(i)/t_01_w_a(i))**(k/(k-1))
	lambda_2w_a(i)=w_2_a(i)/sqrt(2*k*r*t_02_a(i)/(k+1))
	p_02_w_a(i)=p_02_a(i)/(1-(k-1)*lambda_2w_a(i)**2/(k+1))**(k/(k-1))
	omega_sunshi_da(i)=(p_02_w_a(i)-p_02_w_a_ad(i))/(p_11_a(i)-p_01_w_a(i))
	!write(10,*)'静叶出口','i=',i,'t_02_a=',t_02_a(i),'lambda_2_a=',lambda_2_a(i),'p_02_a=',p_02_a(i),'p_03_a=',p_03_a(i),&
	!&'t_03_a=',t_03_a(i),'omega_a=',omega_a(i),'omega_sunshi_ja=',omega_sunshi_ja(i),'lambda_1w_a=',lambda_1w_a(i),&
	!&'t_01_w_a=',t_01_w_a(i),'p_01_w_a=',p_01_w_a(i),'t_02_w_a=',t_02_w_a(i),'p_02_w_a_ad=',p_02_w_a_ad(i),'''lambda_2w_a=',lambda_2w_a(i),&
	!&'p_02_w_a=',p_02_w_a(i),'omega_sunshi_da=',omega_sunshi_da(i) 
    write(*,*)p_02_a(i),p_03_a(i),p_12_a(i)
	end do
	!确定各列叶栅气动计算尺寸
	K_G=0.999
	delta_h=0.5*(1+K_G)
	delta_t=0.5*(1+K_G)
	do i=1,8
	r_1_h(i)=d_h1(i)/2
	r_1_t(i)=d_t1(i)/2
	r_2_h(i)=d_h2(i)/2
	r_2_t(i)=d_t2(i)/2
	r_1_m(i)=(d_h1(i)+d_t1(i))/2
	r_2_m(i)=(d_h2(i)+d_t2(i))/2
	r_he(i)=sqrt(delta_h*r_1_h(i)**2+(1-delta_h)*r_1_t(i)**2)
	r_te(i)=sqrt(delta_t*r_1_t(i)**2+(1-delta_t)*r_1_h(i)**2)
	end do
	!确定径向计算站数目及其位置
	z_r=7
	do i=1,z
	r_r1(1,i)=r_1_h(i)
	do j=1,z_r-1
	r_r1(j+1,i)=r_r1(j,i)+(r_te(i)-r_he(i))/(z_r-1)
	end do
	end do
	!动叶进口截面沿径向气流参数分布规律
	do i=1,z
	A(i)=0.9*2*pi*n*(i-omega_a(i))              !自取
	B(i)=(c_1u_a(i)-A(i)*r_1_m(i))*r_1_m(i)
	  do j=1,z_r
	u_1(j,i)=2*pi*r_r1(j,i)*n/60
	c_1u(j,i)=A(i)*r_r1(1,i)+B(i)/r_r1(1,i)
	c_1a(j,i)=sqrt(c_a1_a(i)**2+2*A(i)**2*(r_1_m(i)**2-r_r1(j,i)**2)+4*A(i)*B(i)*log(r_1_m(i)/r_r1(j,i)))
	c_1(j,i)=sqrt(c_1a(j,i)**2+c_1u(j,i)**2)
	w_1u(j,i)=u_1(j,i)-c_1u(j,i)
	w_1(j,i)=sqrt(c_1a(j,i)**2+w_1u(j,i)**2)
	alpha_1(j,i)=atan(c_1a(j,i)/c_1u(j,i))
	beta_1(j,i)=atan(c_1a(j,i)/w_1u(j,i))
	!计算气体状态（认为径向加工量不变）
	lambda_1(j,i)=c_1(j,i)/sqrt(2*k*r*t_01_a(i)/(k+1))
	t_11(j,i)=t_01_a(i)*(1-(k-1)*lambda_1(j,i)**2/(k+1))
	M_w(j,i)=w_1(j,i)/sqrt(k*r*t_11(j,i))
	sigma_1=(/0.992,0.994,0.996,0.998,0.998,0.997,0.996/)   !校核
	p_01(j,i)=sigma_1(j)*p_01_a(i)
	p_11(j,i)=p_01(j,i)*(1-(k-1)*lambda_1(j,i)**2/(k+1))**(k/(k-1))
	rhoc_11(j,i)=p_11(j,i)/(t_11(j,i)*r)
	end do
	end do
	!流量校核
	G_c_1=0
	do i=1,z
	  do j=1,z_r
	   delta_GA_1(j,i)=rhoc_11(j,i)*c_1a(j,i)
	   if(j>1)then 
	   delta_A_1(j-1,i)=pi*(r_r1(j,i)**2-r_r1(j-1,i)**2)
	   delta_G_1(j-1,i)=delta_A_1(j-1,i)*(delta_GA_1(j,i)-delta_GA_1(j-1,i))/2
	    G_c_1(i)=G_c_1(i)+delta_G_1(j-1,i)
	   end if
	   end do
	  eps_1(i)=(G_c_1(i)-G)/G 
	end do
	do i=1,z
	   t_0zong_1(1)=0
	   p_0zong_1(1)=0
	   do j=1,z_r-1
	   t_0zong_1(i)=delta_G_1(j,i)*(t_01_a(i)+t_01_a(i+1))/2+t_0zong_1(i)
	   p_0zong_1(i)=delta_G_1(j,i)*(p_01_a(i)+p_01_a(i+1))/2+p_0zong_1(i)
	   end do
	   t_0average_1(i)=t_0zong_1(i)/G
	   p_0average_1(i)=p_0zong_1(i)/G
	end do 
	 !动叶出口截面
	 do i=1,z
	   E(i)=(1-(d_h1(i)-d_t1(i))**2)*d_t1(i)**2/((1-(d_h2(i)-d_t2(i))**2)*d_t2(i)**2)
	    do j=1,z_r
		r_r2(j,i)=r_2_m(i)*sqrt(1+(r_1_m(i)/r_1_m(i))**2*((r_r1(j,i)/r_1_m(i))**2-1)/E(i))
		u_2(j,i)=2*pi*r_r2(j,i)*n/60
		h(i)=u_2_a(i)*c_2u_a(i)-u_2_a(i)*c_2u_a(i)
		C(i)=E(i)*A(i)
		D(i)=h(i)/(2*pi*n)+A(i)*r_1_m(i)**2-c_2_a(i)**2+B(i)
		c_2u(j,i)=C(i)*r_r2(j,i)+D(i)/r_r2(j,i)
		c_2a(j,i)=sqrt(c_a2_a(i)**2+2*C(i)**2*(r_2_m(i)**2-r_r2(j,i)**2)+4*C(i)*D(i)*log(r_2_m(i)/r_r2(j,i)))
		c_2(j,i)=sqrt(c_2a(j,i)**2+c_2u(j,i)**2)
		w_2u(j,i)=u_2(j,i)-c_2u(j,i)
		w_2(j,i)=sqrt(c_2a(j,i)**2+w_2u(j,i)**2)
		alpha_2(j,i)=atan(c_2a(j,i)/c_2u(j,i))
	    beta_2(j,i)=atan(c_2a(j,i)/w_2u(j,i))
	    end do
	end do
	!计算气体状态（动叶）
	  eps_2=(/0.996,0.996,0.995,0.995,0.995,0.993,0.990/)
	do i=1,z
	delta_t0(i)=h(i)/c_p
	t_02(i)=t_01_a(i)+delta_t0(i)
	  do j=1,z_r
	  lambda_2(j,i)=c_2(j,i)/sqrt(2*k*r*t_02(i)/(k+1))
	  t_12(j,i)=t_02(i)*(1-(k-1)*lambda_2(j,i)**2/(k+1))
	  M_c(j,i)=c_2(j,i)/sqrt(k*r*t_12(j,i))
	  t_01_w(j,i)=t_11(j,i)+w_1(j,i)**2/(2*c_p)
	  lambda_1w(j,i)=w_1(j,i)/sqrt(2*k*r*t_01_w(j,i)/(k+1))
	  p_01_w(j,i)=p_11(j,i)/(1-(k-1)*lambda_1w(j,i)**2/(k+1))**(k/(k-1))
	  t_02_w(j,i)=t_12(j,i)+w_2(j,i)**2/(2*c_p)
	  lambda_2w(j,i)=w_2(j,i)/sqrt(2*k*r*t_02_w(j,i)/(k+1))
	  p_02_w_ad(j,i)=p_01_w(j,i)*(t_02_w(j,i)/t_01_w(j,i))**(k*(k-1))
	  p_02_w(j,i)= p_02_w_ad(j,i)*eps_2(j)
	  p_12(j,i)=p_02_w(j,i)*(1-(k-1)*lambda_2w(j,i)**2/(k+1))**(k/(k-1))
	  p_02(j,i)=p_12(j,i)/(1-(k-1)*lambda_2(j,i)**2/(k+1))**(k/(k-1))
	  rhoc_12(j,i)=p_12(j,i)/(t_12(j,i)*r)
	  delta_GA_2(j,i)=rhoc_12(j,i)*c_2a(j,i)
	  end do 
	end do
	G_c_2=0
	do i=1,z
	  do j=1,z_r-1
	  delta_A_2(j,i)=pi*(r_r2(j+1,i)**2-r_r2(j,i)**2)
	  delta_G_2(j,i)=delta_A_2(j,i)*(delta_GA_2(j+1,i)+delta_GA_2(j,i))/2
	  G_c_2(i)=G_c_2(i)+delta_G_2(j,i)
	  end do
	  eps_2(i)=(G_c_2(i)-G)/G
	end do 
	do i=1,z
       do j=1,z_r
	   omega(j,i)=1-(c_1u(j,i)+c_2u(j,i))/(u_1(j,i)+u_2(j,i))
	   eta_adi(j,i)=((p_01(j,i+1)/p_02(j,i))**((k-1)/k)-1)/(t_02(i)/t_01_a(i)-1)
	   end do
	end do
	do i=1,z
	   t_0zong_2(1)=0
	   p_0zong_2(1)=0
	   eta_zong(1)=0
	   do j=1,z_r-1
	   t_0zong_2(i)=delta_G_2(j,i)*(t_02(i)+t_02(i))/2+t_0zong_2(i)
	   p_0zong_2(i)=delta_G_2(j,i)*(p_02(j+1,i)+p_02(j,i))/2+p_0zong_2(i)
	   eta_zong(i)=delta_G_2(j,i)*(eta_adi(j,i)+eta_adi(j+1,i))/2+eta_zong(i)
	   end do
	   t_0average_2(i)=t_0zong_2(i)/G
	   p_0average_2(i)=p_0zong_2(i)/G
	  eta_average(i)=eta_zong(i)/G
	end do 
	close(10)
    stop
    end program  one
	
	!#####################截面参数子程序###################################
	subroutine  jmcs(c1_a1,d1_1,alpha1_1,t1_0,p1_0,p1_1,t1_1,rho1_1,c1_1,q1_1,lambda1_1,d1_t1,d1_h1,a1_1,L1_1)
	real::c1_a1,d1_1,alpha1_1,t1_0,p1_0                !自取
	real::p1_1,t1_1,rho1_1,c1_1,q1_1,lambda1_1,d1_t1,d1_h1,a1_1,L1_1  !求解
	c1_1=c1_a1/sin(alpha1_1*pi/180)
	lambda1_1=c1_1/sqrt(2*r*k*t1_0/(k+1))
	t1_1=(1-(k-1)*lambda1_1**2/(k+1))*t1_0
	p1_1=(1-(k-1)*lambda1_1**2/(k+1))**(k/(k-1))*p1_0
	q1_1=lambda1_1*((k+1)/2-(k-1)*(lambda1_1**2)/2)**(1/(k-1))
	a1_1=G*sqrt(t1_0)/(k1*q1_1*p1_0*sin(alpha1_1*pi/180))
	d1_t1=sqrt(4*a1_1/(pi*(1-d1_1**2)))
	d1_h1=d1_t1*d1_1
	L1_1=(d1_t1-d1_h1)/2
	return
	end subroutine jmcs
	
	
	
	
	
	
	
	
	
	
	