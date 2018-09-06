classdef MieScattering
       
       properties
            value_N
            value_radius
            value_Eo
            value_ko
            value_ns
       end
       
       methods
                function obj = MieScattering(inp_N,inp_radius,inp_Eo,inp_ko,inp_ns)
                    if (isscalar(inp_N)==0)
                        error('Error: inp_N must be a scalar');
                    end
                    if (isscalar(inp_radius)==0)
                        error('Error: inp_radius must be a scalar');
                    end
                    if (isscalar(inp_Eo)==0)
                        error('Error: inp_Eo must be a scalar');
                    end
                    if (isscalar(inp_ko)==0)
                        error('Error: inp_ko must be a scalar');
                    end
                    if (isscalar(inp_ns)==0)
                        error('Error: inp_ns must be a scalar');
                    end
                    obj.value_N = inp_N;
                    obj.value_radius = inp_radius;
                    obj.value_Eo = inp_Eo;
                    obj.value_ko = inp_ko;
                    obj.value_ns = inp_ns;
                end
                
                function[sbi] = BistaticRCS(obj,phi,theta)
                    if (isscalar(phi)==0)
                        error('Error: phi must be a scalar');
                    end
                    N = obj.value_N;
                    radius = obj.value_radius;
                    ko = obj.value_ko;
                    ns = obj.value_ns;
                    z = radius*ko;
                    ftheta = zeros(1,length(theta));
                    fphi = zeros(1,length(theta));
                    for n=1:N
                        if abs(ns) == 0
                            [A, B] = MieScattering.Parameters_PEC_Scattered(n,z);
                        else
                            [A, B] = MieScattering.Parameters_Dielectric_Scattered(n,z,ns);
                        end
                        aLfn = MieScattering.Associated_Legendre_Functions(n,theta);
                        aLf = aLfn(2,:)./sin(theta);
                        daLf = MieScattering.Derrivative_Associated_Legendre_Functions(1,n,theta);
                        i = find(theta==0);
                        aLf(i) = -(n*(n+1))/2;
                        daLf(i) = aLf(i);
                        i = find(theta==pi);
                        aLf(i) = (-1)^n*(n*(n+1))/2;
                        daLf(i) = -aLf(i);
                        i = find(theta==-pi);
                        aLf(i) = (-1)^n*(n*(n+1))/2;
                        daLf(i) = -aLf(i);
                        ftheta = ftheta + (2*n+1)/(n*(n+1))*(A.*aLf + B.*daLf);
                        fphi = fphi + (2*n+1)/(n*(n+1))*(A.*daLf + B.*aLf);
                    end
                    ftheta = -1i*1/ko*cos(phi)*ftheta;
                    fphi = 1i*1/ko*sin(phi)*fphi;
                    sd = abs(ftheta).^2 + abs(fphi).^2;
                    sbi = 4*pi*sd;
                end
              
                function[smo] = MonostaticRCS(obj)
                    N = obj.value_N;
                    radius = obj.value_radius;
                    ko = obj.value_ko;
                    ns = obj.value_ns;
                    z = radius*ko;
                    smo = zeros(1,length(z));
                    for n=1:N
                        if abs(ns) == 0
                            [A, B] = MieScattering.Parameters_PEC_Scattered(n,z);
                        else
                            [A, B] = MieScattering.Parameters_Dielectric_Scattered(n,z,ns);
                        end
                        smo = smo + (-1)^n*(2*n+1)*(A-B);
                    end
                    smo = 1./(z.^2).*abs(smo).^2;
                end 
              
                function[EincTotal, HincTotal] = Incident_Fields(obj,phi,theta,observation)
                    if (isscalar(phi)==0 || isscalar(observation)==0)
                        error('Error: phi and observation point must be scalars');
                    end
                    N = obj.value_N;
                    Eo = obj.value_Eo;
                    ko = obj.value_ko;
                    c = 299792458;
                    mu = 4*pi*10^(-7);
                    epsilon = 1/c^2/mu;
                    r = observation*ko;   
                    Zo = sqrt(mu/epsilon);
                    Einctheta = zeros(length(theta),1)';
                    Eincphi = zeros(length(theta),1)';
                    Eincrho = zeros(length(theta),1)';
                    Hinctheta = zeros(length(theta),1)';
                    Hincphi = zeros(length(theta),1)';
                    Hincrho = zeros(length(theta),1)';
                    for n=1:N
                        [me, mo] = MieScattering.Spherical_Wave_function_m(1,r,theta,phi,1,n);
                        [ne, no] = MieScattering.Spherical_Wave_function_n(1,r,theta,phi,1,n);
                        Einctheta = Einctheta + (-1i)^n*(2*n+1)/(n*(n+1))*(mo(1,:) + 1i*ne(1,:));
                        Hinctheta = Hinctheta + (-1i)^n*(2*n+1)/(n*(n+1))*(me(1,:) - 1i*no(1,:));
                        Eincphi = Eincphi + (-1i)^n*(2*n+1)/(n*(n+1))*(mo(2,:) + 1i*ne(2,:));
                        Hincphi = Hincphi + (-1i)^n*(2*n+1)/(n*(n+1))*(me(2,:) - 1i*no(2,:));
                        Eincrho = Eincrho + (-1i)^n*(2*n+1)/(n*(n+1))*(mo(3,:) + 1i*ne(3,:));
                        Hincrho = Hincrho + (-1i)^n*(2*n+1)/(n*(n+1))*(me(3,:) - 1i*no(3,:));
                    end
                    EincTotal = -Eo*[Eincrho; Einctheta; Eincphi];
                    HincTotal = Eo/Zo*[Hincrho; Hinctheta; Hincphi];
                end
              
                function[EscaTotal, HscaTotal] = Interior_Fields_Line(obj,x,y,z)
                    if isscalar(x) == 0
                        len = length(x);
                    elseif isscalar(y) == 0 
                        len = length(y);
                    elseif isscalar(z) == 0
                        len = length(z);
                    end
                    if ~((len/length(x) == 1 || len/length(x) == len) && (len/length(y) == 1 || len/length(y) == len) && (len/length(z) == 1 || len/length(z) == len))
                        error('Error: x,y,z must have the same length or be scalars');
                    end
                    [phi,theta,observation] = cart2sph(x,y,z);
                    theta = theta-pi/2;
                    if isscalar(theta) == 1
                        theta = theta*ones(len,1);
                    end
                    if isscalar(phi) == 1 
                        phi = phi*ones(len,1);
                    end
                    if isscalar(observation) == 1 
                        observation = observation*ones(len,1);
                    end
                    N = obj.value_N;
                    radius = obj.value_radius;
                    Eo = obj.value_Eo;
                    ko = obj.value_ko;
                    ns = obj.value_ns;
                    c = 299792458;
                    mu = 4*pi*10^(-7);
                    epsilon = 1/c^2/mu;
                    r = observation*ko;
                    z = radius*ko;
                    Zo = sqrt(mu/epsilon);
                    Ex = zeros(len,1);
                    Ey = zeros(len,1);
                    Ez = zeros(len,1);
                    Hx = zeros(len,1);
                    Hy = zeros(len,1);
                    Hz = zeros(len,1);
                    for i=1:len
                        Escatheta = 0;
                        Escaphi = 0;
                        Escarho = 0;
                        Hscatheta = 0;
                        Hscaphi = 0;
                        Hscarho = 0;
                        for n=1:N
                            [C, D] = MieScattering.Parameters_Dielectric_Interior(n,z,ns);
                            [me, mo] = MieScattering.Spherical_Wave_function_m(1,ns*r(i),theta(i),phi(i),1,n);
                            [ne, no] = MieScattering.Spherical_Wave_function_n(1,ns*r(i),theta(i),phi(i),1,n);
                            Escatheta = Escatheta + (-1i)^n*(2*n+1)/(n*(n+1))*(C*mo(1,:) + 1i*D*ne(1,:));
                            Hscatheta = Hscatheta + (-1i)^n*(2*n+1)/(n*(n+1))*(D*me(1,:) - 1i*C*no(1,:));
                            Escaphi = Escaphi + (-1i)^n*(2*n+1)/(n*(n+1))*(C*mo(2,:) + 1i*D*ne(2,:));
                            Hscaphi = Hscaphi + (-1i)^n*(2*n+1)/(n*(n+1))*(D*me(2,:) - 1i*C*no(2,:));
                            Escarho = Escarho + (-1i)^n*(2*n+1)/(n*(n+1))*(C*mo(3,:) + 1i*D*ne(3,:));
                            Hscarho = Hscarho + (-1i)^n*(2*n+1)/(n*(n+1))*(D*me(3,:) - 1i*C*no(3,:));
                        end
                        Ex(i) = sin(theta(i))*cos(phi(i))*Escarho + cos(theta(i))*cos(phi(i))*Escatheta - sin(phi(i))*Escaphi;
                        Ey(i) = sin(theta(i))*sin(phi(i))*Escarho + cos(theta(i))*sin(phi(i)).*Escatheta + cos(phi(i))*Escaphi;
                        Ez(i) = cos(theta(i)).*Escarho - sin(theta(i))*Escatheta;
                        Hx(i) = sin(theta(i))*cos(phi(i)).*Hscarho + cos(theta(i))*cos(phi(i))*Hscatheta - sin(phi(i)).*Hscaphi;
                        Hy(i) = sin(theta(i))*sin(phi(i)).*Hscarho + cos(theta(i))*sin(phi(i))*Hscatheta + cos(phi(i)).*Hscaphi;
                        Hz(i) = cos(theta(i))*Hscarho - sin(theta(i)).*Hscatheta;
                    end
                    EscaTotal(:,1) = Eo*Ex;
                    EscaTotal(:,2) = Eo*Ey;
                    EscaTotal(:,3) = Eo*Ez;
                    HscaTotal(:,1) = Eo*ns/Zo*Hx;
                    HscaTotal(:,2) = Eo*ns/Zo*Hy;
                    HscaTotal(:,3) = Eo*ns/Zo*Hz;
                end
              
                function[EscaTotal,HscaTotal] = Scattered_Fields(obj,phi,theta,observation)
                    if (isscalar(phi)==0 || isscalar(observation)==0)
                        error('Error: phi and observation point must be scalars');
                    end
                    N = obj.value_N;
                    radius = obj.value_radius;
                    Eo = obj.value_Eo;
                    ko = obj.value_ko;
                    ns = obj.value_ns;
                    c = 299792458;
                    mu = 4*pi*10^(-7);
                    epsilon = 1/c^2/mu;
                    r = observation*ko;
                    z = radius*ko;
                    Zo = sqrt(mu/epsilon);
                    Escatheta = zeros(length(theta),1)';
                    Escaphi = zeros(length(theta),1)';
                    Escarho = zeros(length(theta),1)';
                    Hscatheta = zeros(length(theta),1)';
                    Hscaphi = zeros(length(theta),1)';
                    Hscarho = zeros(length(theta),1)';
                    for n=1:N
                        if abs(ns) == 0
                            [A, B] = MieScattering.Parameters_PEC_Scattered(n,z);
                        else
                            [A, B] = MieScattering.Parameters_Dielectric_Scattered(n,z,ns);
                        end
                        [me, mo] = MieScattering.Spherical_Wave_function_m(4,r,theta,phi,1,n);
                        [ne, no] = MieScattering.Spherical_Wave_function_n(4,r,theta,phi,1,n);
                        Escatheta = Escatheta + (-1i)^n*(2*n+1)/(n*(n+1))*(A*mo(1,:) + 1i*B*ne(1,:));
                        Hscatheta = Hscatheta + (-1i)^n*(2*n+1)/(n*(n+1))*(B*me(1,:) - 1i*A*no(1,:));
                        Escaphi = Escaphi + (-1i)^n*(2*n+1)/(n*(n+1))*(A*mo(2,:) + 1i*B*ne(2,:));
                        Hscaphi = Hscaphi + (-1i)^n*(2*n+1)/(n*(n+1))*(B*me(2,:) - 1i*A*no(2,:));
                        Escarho = Escarho + (-1i)^n*(2*n+1)/(n*(n+1))*(A*mo(3,:) + 1i*B*ne(3,:));
                        Hscarho = Hscarho + (-1i)^n*(2*n+1)/(n*(n+1))*(B*me(3,:) - 1i*A*no(3,:));
                    end
                    EscaTotal = -Eo*[Escarho; Escatheta; Escaphi];
                    HscaTotal = Eo/Zo*[Hscarho; Hscatheta; Hscaphi];
                end
              
                function[Q_scat,Q_ext] = Power_Efficiency(obj)
                    N = obj.value_N;
                    radius = obj.value_radius;
                    Eo = obj.value_Eo;
                    ko = obj.value_ko;
                    ns = obj.value_ns;
                    z = radius*ko;
                    G = pi*radius^2;
                    C_ext = zeros(1,length(z));
                    C_scat = zeros(1,length(z));
                    for n=1:N
                        [A, B] = MieScattering.Parameters_Dielectric_Scattered(n,z,ns);
                        C_ext = C_ext + (2*n+1)*(A+B);
                        C_scat = C_scat + (2*n+1)*(abs(A).^2+abs(B).^2);
                    end
                    C_ext = 2*pi/ko^2*real(C_ext);
                    Q_ext = C_ext./G/Eo;
                    C_scat = 2*pi/ko^2*C_scat;
                    Q_scat = C_scat./G/Eo;
                end
       end
       
       methods (Static)
                function[aLf] = Associated_Legendre_Functions(n,theta)
                    aLf = legendre(n,cos(theta));
                end            
                function[daLf] = Derrivative_Associated_Legendre_Functions(m,n,theta)
                    aLf = legendre(n,cos(theta));
                    aLfn1 = legendre(n+1,cos(theta));
                    daLf = ((n-m+1).*aLfn1(m+1,:) - (n+1)*cos(theta).*aLf(m+1,:))./(sin(theta));
                end
                function[Sb1k] = Spherical_Bessel_function_1k(n,z)
                    if (z==0 && n==0)
                        Sb1k = 1;
                    elseif (z==0 && n~=0)
                        Sb1k = 0;
                    else
                        Sb1k = sqrt(pi./(2*z)).*besselj(n+1/2,z);
                    end
                end
                function[Sb2k] = Spherical_Bessel_function_2k(n,Z)
                    Sb2k = sqrt(pi./(2*Z)).*bessely(n+1/2,Z);
                end
                function[Shnk] = Spherical_Hankel_function(n,Z,k)
                    Shnk = sqrt(pi./(2*Z)).*besselh(n+1/2,k,Z);
                end
                function[dSb] = Differentation_Spherical_Bessel_function(n,Z,k)
                    if k==1
                        zk = MieScattering.Spherical_Bessel_function_1k(n,Z);
                        zkn1 = MieScattering.Spherical_Bessel_function_1k(n-1,Z);
                    elseif k==2
                        zk = MieScattering.Spherical_Bessel_function_2k(n,Z);
                        zkn1 = MieScattering.Spherical_Bessel_function_2k(n-1,Z);
                    elseif k==3
                        zk = MieScattering.Spherical_Hankel_function(n,Z,1);
                        zkn1 = MieScattering.Spherical_Hankel_function(n-1,Z,1);
                    elseif k==4
                        zk = MieScattering.Spherical_Hankel_function(n,Z,2);
                        zkn1 = MieScattering.Spherical_Hankel_function(n-1,Z,2);
                    else
                        error('Wrong kind of Spherical Bessel or Hankel function') 
                    end
                    dSb = Z.*zkn1 - n*zk;
                end
                function[me, mo] = Spherical_Wave_function_m(k,z,theta,phi,m,n)
                    if k==1
                        zk = MieScattering.Spherical_Bessel_function_1k(n,z);
                    elseif k==2
                        zk = MieScattering.Spherical_Bessel_function_2k(n,z);
                    elseif k==3
                        zk = MieScattering.Spherical_Hankel_function(n,z,1);
                    elseif k==4
                        zk = MieScattering.Spherical_Hankel_function(n,z,2);
                    else
                        error('Wrong kind of Spherical Bessel or Hankel function') 
                    end
                    aLfn = MieScattering.Associated_Legendre_Functions(n,theta);
                    aLf = aLfn(m+1,:)./sin(theta);
                    daLf = MieScattering.Derrivative_Associated_Legendre_Functions(m,n,theta);
                    i = find(theta==0);
                    aLf(i) = -(n*(n+1))/2;
                    daLf(i) = aLf(i);
                    i = find(theta==pi);
                    aLf(i) = (-1)^n*(n*(n+1))/2;
                    daLf(i) = -aLf(i);
                    i = find(theta==-pi);
                    aLf(i) = (-1)^n*(n*(n+1))/2;
                    daLf(i) = -aLf(i);
                    me(1,:) = -zk.*aLf*m.*sin(m.*phi);
                    me(2,:) = -zk.*daLf.*cos(m.*phi);
                    me(3,:) = zeros(length(theta),1);
                    mo(1,:) = zk.*aLf.*m.*cos(m.*phi);
                    mo(2,:) = -zk.*daLf.*sin(m.*phi);
                    mo(3,:) = zeros(length(theta),1);
                end
                function[ne, no] = Spherical_Wave_function_n(k,z,theta,phi,m,n)
                    if k==1
                        zkz = (MieScattering.Spherical_Bessel_function_1k(n-1,z) + MieScattering.Spherical_Bessel_function_1k(n+1,z))/(2*n+1);
                        dsBz = MieScattering.Spherical_Bessel_function_1k(n-1,z) - n/(2*n+1)*(MieScattering.Spherical_Bessel_function_1k(n-1,z) + MieScattering.Spherical_Bessel_function_1k(n+1,z));
                    elseif k==2
                        zkz = (MieScattering.Spherical_Bessel_function_2k(n-1,z) + MieScattering.Spherical_Bessel_function_2k(n+1,z))/(2*n+1);
                        dsBz = MieScattering.Spherical_Bessel_function_2k(n-1,z) - n/(2*n+1)*(MieScattering.Spherical_Bessel_function_2k(n-1,z) + MieScattering.Spherical_Bessel_function_2k(n+1,z));
                    elseif k==3
                        zkz = (MieScattering.Spherical_Hankel_function(n-1,z,1) + MieScattering.Spherical_Hankel_function(n+1,z,1))/(2*n+1);
                        dsBz = MieScattering.Spherical_Hankel_function(n-1,z,1) - n/(2*n+1)*(MieScattering.Spherical_Hankel_function(n-1,z,1) + MieScattering.Spherical_Hankel_function(n+1,z,1));
                    elseif k==4
                        zkz = (MieScattering.Spherical_Hankel_function(n-1,z,2) + MieScattering.Spherical_Hankel_function(n+1,z,2))/(2*n+1);
                        dsBz = MieScattering.Spherical_Hankel_function(n-1,z,2) - n/(2*n+1)*(MieScattering.Spherical_Hankel_function(n-1,z,2) + MieScattering.Spherical_Hankel_function(n+1,z,2));
                    else
                        error('Wrong kind of Spherical Bessel or Hankel function') 
                    end
                    aLfn = MieScattering.Associated_Legendre_Functions(n,theta);
                    aLf = aLfn(m+1,:)./sin(theta);
                    daLf = MieScattering.Derrivative_Associated_Legendre_Functions(m,n,theta);
                    i = find(theta==0);
                    aLf(i) = -(n*(n+1))/2;
                    daLf(i) = aLf(i);
                    i = find(theta==pi);
                    aLf(i) = (-1)^n*(n*(n+1))/2;
                    daLf(i) = -aLf(i);
                    i = find(theta==-pi);
                    aLf(i) = (-1)^n*(n*(n+1))/2;
                    daLf(i) = -aLf(i);
                    ne(1,:) = dsBz.*daLf.*cos(m.*phi);
                    ne(2,:) = -dsBz.*aLf.*m.*sin(m.*phi);
                    ne(3,:) = n.*(n+1).*zkz.*aLf.*sin(theta).*cos(m*phi);
                    no(1,:) = dsBz.*daLf.*sin(m.*phi);
                    no(2,:) = dsBz.*aLf.*m.*cos(m.*phi);
                    no(3,:) = n.*(n+1).*zkz.*aLf.*sin(theta).*sin(m.*phi);
                end
                function[A, B] = Parameters_Dielectric_Scattered(n,z,ns)
                    A = -(MieScattering.Spherical_Bessel_function_1k(n,z*ns).*MieScattering.Differentation_Spherical_Bessel_function(n,z,1) - MieScattering.Spherical_Bessel_function_1k(n,z).*MieScattering.Differentation_Spherical_Bessel_function(n,z*ns,1))./(MieScattering.Spherical_Bessel_function_1k(n,z*ns).*MieScattering.Differentation_Spherical_Bessel_function(n,z,4) - MieScattering.Spherical_Hankel_function(n,z,2).*MieScattering.Differentation_Spherical_Bessel_function(n,z*ns,1));
                    B = -(MieScattering.Spherical_Bessel_function_1k(n,z).*MieScattering.Differentation_Spherical_Bessel_function(n,z*ns,1) - ns^2*MieScattering.Spherical_Bessel_function_1k(n,z*ns).*MieScattering.Differentation_Spherical_Bessel_function(n,z,1))./(MieScattering.Spherical_Hankel_function(n,z,2).*MieScattering.Differentation_Spherical_Bessel_function(n,z*ns,1) - ns^2*MieScattering.Spherical_Bessel_function_1k(n,z*ns).*MieScattering.Differentation_Spherical_Bessel_function(n,z,4));
                end
                function[C, D] = Parameters_Dielectric_Interior(n,z,ns)
                    C = -(-MieScattering.Spherical_Bessel_function_1k(n,z).*MieScattering.Differentation_Spherical_Bessel_function(n,z,4) + MieScattering.Spherical_Hankel_function(n,z,2).*MieScattering.Differentation_Spherical_Bessel_function(n,z,1))./(MieScattering.Spherical_Bessel_function_1k(n,z*ns).*MieScattering.Differentation_Spherical_Bessel_function(n,z,4) - MieScattering.Spherical_Hankel_function(n,z,2).*MieScattering.Differentation_Spherical_Bessel_function(n,z*ns,1));
                    D = -ns*(MieScattering.Spherical_Bessel_function_1k(n,z).*MieScattering.Differentation_Spherical_Bessel_function(n,z,4) - MieScattering.Spherical_Hankel_function(n,z,2).*MieScattering.Differentation_Spherical_Bessel_function(n,z,1))./(MieScattering.Spherical_Hankel_function(n,z,2).*MieScattering.Differentation_Spherical_Bessel_function(n,z*ns,1) - ns^2*MieScattering.Spherical_Bessel_function_1k(n,z*ns).*MieScattering.Differentation_Spherical_Bessel_function(n,z,4));
                end
                function[A, B] = Parameters_PEC_Scattered(n,z)
                    A = -(MieScattering.Spherical_Bessel_function_1k(n,z))./(MieScattering.Spherical_Hankel_function(n,z,2));
                    B = -(MieScattering.Differentation_Spherical_Bessel_function(n,z,1))./(MieScattering.Differentation_Spherical_Bessel_function(n,z,4));
                end        
       end  
       
end
