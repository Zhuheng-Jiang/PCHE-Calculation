%=====================================%
% This program is designed to caculate the heat transfer performance and the pressure drop of the hybrid PCHE, for the FLNG technology.
% And this code was produced by Zhuheng Jiang for his own Bachelor's graduation project in Xi'an Jiaotong University.

%The following project is based on the MATLAB

%Define size in meter 
diahx=1.5e-3;
thihx=1.35e-3;
dfin=5e-3;
tfin=1e-4;
wide=2e-3;
channels=150;
L=0.4218;  %initial length in meters; threshold
finnum=1;
% downfg=0;%set two bool to avoid endless loop
% upfg=0;
% for Iter=1:1:50
    
    %This part defines needed parameters 
    
    parts=60;
    lunit=L/parts; %Length of every calculation unit(divided into 50 uniform part)
    lcharhx=4*diahx*pi/(pi+2); %Define the characteristic length
%     lcharfin= 4*(((wide-finnum*tfin)/(finnum+1))*(dfin-tfin))/((dfin-tfin)*2+(wide-finnum*tfin)*2/finnum); % Same
     lcharfin= 4*(((wide-tfin))*(dfin-tfin)/2)/((dfin-tfin)*2+(wide-tfin)); % Same
    %=============================================================%
    %count the area(Both sides)
    Ahx=lunit*pi*diahx/2+diahx*lunit;
    Ahxcross=pi*diahx*diahx/8;               %cross section area

    Afin=lunit*(2*(wide-finnum*tfin)+finnum*(dfin-tfin)*2);  %Efficiency of fins
    Afincross=(wide-finnum*tfin)*(dfin-tfin);

    Asteel=lunit*wide;         %====use the wide
    thisteel=((thihx-2e-4)*diahx-diahx^2*pi/8)/diahx;   %Average
    ksteel=16.3;    %Heat conductivity of 316L in W/(m*K) 
    %=====================================================%

    %Give the inlet(hot) and outlet(cold) temp as initial condition
    Tin=300.0;
    tout=171.4;   %outlet's approximate temp of the cold path
    Tpse=215.16;  %Pseudo-Critical
    Tcri=200.01;

    PNG=8000;      %unit in kPa(REFPROP's pressure unit is kPa)   8000kPa
    PNit=500;     %0.5MPa
    mNG=0.012/channels;    %Flow rate in kg/s
    mNit=0.118/channels;

    Tu=Tin;       %The initial temp of the u=1 unit
    tu=tout;
    %=======define three 1*51 arrays to store the Tout, tout & q
    Tarray= zeros(1,parts+1);
    tarray= zeros(1,parts+1);
    qarray=zeros(1,parts);
    Q=0;

    Tarray(1)=Tin;
    tarray(1)=tout;

    PNGtotal=0;
    PNittotal=0;
    %Write the loop
    %ReNitarray=zeros(parts);
    for unit = 1 : 1: parts;

        %Twall= ！！！！！！！！！！！！！！

        %======================================================================
        Twall=(Tu+tu)/2;   %In every unit the initial Twall
        %==========1. Heat condition==============================
        %====Consult the REFPROP===
        kNG=refpropm('L', 'T', Tu, 'P', PNG, 'NGSAMPLE.mix' );   % W/k*m
        kNit=refpropm('L', 'T', tu, 'P', PNit, 'NITROGEN' );
        CpNG=refpropm('C','T', Tu, 'P', PNG, 'NGSAMPLE.mix' );   %  J/kg*K
        CpNit=refpropm('C', 'T', tu, 'P', PNit, 'NITROGEN' );
        roNG=refpropm('D', 'T', Tu, 'P', PNG, 'NGSAMPLE.mix');   %Density 
        roNit=refpropm('D', 'T', tu, 'P', PNit, 'NITROGEN');
        miuNG=refpropm('V', 'T', Tu, 'P', PNG, 'NGSAMPLE.mix');  %Dynamic viscosity
        miuNit=refpropm('V', 'T', tu, 'P', PNit, 'NITROGEN');
        PrNG=refpropm('^', 'T', Tu, 'P', PNG, 'NGSAMPLE.mix' );  %Prantl number
        PrNit=refpropm('^', 'T', tu, 'P', PNit, 'NITROGEN');
        ReNG= mNG*lcharhx/(miuNG*Ahxcross);     %Viscosity will change
        ReNit=mNit*lcharfin/(miuNit*Afincross);
        %ReNitarray(unit)=ReNit;

        for iter=1:1:50
            %======Calculate the NuNG for PCHE side===========
            miuNGwall=refpropm('V', 'T', Twall, 'P', PNG, 'NGSAMPLE.mix');  
            roNGwall=refpropm('D', 'T', Twall, 'P', PNG, 'NGSAMPLE.mix');
            %Different section has different Nu
            if Twall>=Tpse       %Supercritical section    
                a=0.0073;
                b=0.9463;
                c=0.8037;
                d=-0.2814;
            elseif Twall>=Tcri   %Pseudo-critical section
                a=0.0069;
                b=0.9619;
                c=0.2613;
                d=-0.9011;
            else              %Subcritical section
                a=0.00054;
                b=1.3384;
                c=-0.0585;
                d=0;
            end

            NuNG= a* (ReNG^b)*(PrNG^c)*(roNG/roNGwall)^d;    %roWall!!！！！！！！！！！！

            %====NuNit=====  
            ct=(tu/Twall)^0.45;
            f=(1.82*log10(ReNit)-1.64)^-2;      %Darcy Friction
            NuNit=(((f/8)*(ReNit-1000)*PrNit)/(1+12.7*(f/8)^0.5*(PrNit^1.5-1)))*(1+(lcharfin/lunit)^1.5)*ct;    %NItrogen is heated, use the D-B 
            %===The wavy plain fin has to count the efficiency
            %ei=()
            %Effifin=

            %====Heat convetion coefficient===
            hNG=NuNG* kNG/ lcharhx;  %NG
            hNit=NuNit* kNit/ lcharfin;  %N2
            Twallnew=((hNG*Tu*Ahx)+(hNit*tu*Afin))/((hNG*Ahx)+(hNit*Afin)); %Calculate the new Twall

            if abs(Twallnew-Twall) > 1e-3   %Twall and Twallnew are the same
                Twall=Twallnew;
            else
                break;
            end
        end

        %====Overall Heat resistance===
%         effm=((2*hNit/(ksteel*tfin))*(1+tfin/(wide/finnum)))^0.5;
%         effb=dfin/2;
%         efffin=tanh(effm*effb)/(effm*effb);
%         
        efffin=1;

        Rt = 1/(hNG*Ahx)+ 1/(efffin*hNit*Afin)+ thisteel/(ksteel * Asteel);
        qu= 2*(Tu-tu) / (2*Rt +1/(mNit*CpNit)-1/(mNG*CpNG));
        Tout= Tu- qu/(mNG*CpNG) ;  
        tout= tu- qu/(mNit*CpNit); 

        %====Renew the value to the new Tu and tout====

        Tu=Tout;
        tu=tout;
        Tarray(unit+1)=Tout;
        tarray(unit+1)=tout;
        qarray(unit)=qu;
        Q=Q+qu;

        %================================================================
        %=============2. Pressure drop======================
        % 1. NG side PCHE pressure reduction
        if Twall>=Tpse       %Supercritical section
            af=0.9041;
            bf=-0.2747;
            fiso=0.25*(1.82*log10(ReNG)-1.64)^(-2);
            fNG=a*fiso*(miuNGwall/miuNG)^b;
        elseif Twall>=Tcri   %Pseudo-critical section
            af=1;
            bf=-0.3132;
            fiso=0.25*(1.82*log10(ReNG)-1.64)^(-2);
            fNG=a*fiso*(miuNGwall/miuNG)^b;

        else              %Subcritical section
            af=4.2717;
            bf=-0.7828;
            cf=0.3098;
            fNG=a*ReNG^b*PrNG^c;
        end
        PreduNG = 2*(lunit/lcharhx)*(mNG^2/(roNG*Ahxcross^2))*fNG;   %Fanning here

        PNGtotal=PNGtotal + PreduNG;
        %2. Nitrogen pressure reduction
        fNit= (1.82*log10(ReNit)-1.64)^-2;
        PreduNit = 0.5*(lunit/lcharfin)*(mNit^2/(roNit*Afincross^2))*fNit;    %Darcy here
        PNittotal=PNittotal + PreduNit;
    end
    
    %==3. This part is the judgment for deciding the best Length in 
    %(Now this part is futile, all the check procedures need to be done manually)
    %=======given channels=======
%     if abs(Tarray(51)-110)>1
%         if Tarray(51)<110    
%             if downfg==1&&upfg==1
%                 Lfinal=L-0.1;   %本应该减去0.2，由于已经进行过一次循环，直接取中间值
%                 break;
%             else
%                 downfg=1;
%                 L=L-0.2;
%             end
%         elseif Tarray(51)>110
%             if downfg==1&&upfg==1
%                 Lfinal=L+0.1;   %本应该加上0.2，由于已经进行过一次循环，直接取中间值
%                 break;
%             else
%                 upfg=1;
%                 L=L+0.2;
%             end
%         end
%     else
%         Lfinal=L;
%         break;
%     end
            
    
% end




