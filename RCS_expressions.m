%%Radiofrequency Circuits and Systems useful expressions
function RCS_expressions(option)
    remain= 1;
    choice= 0;
    while remain, 
    if nargin() > 1, 
        disp("Invalid number of arguments\n"),
    elseif nargin == 1 || choice == 1,
        if choice == 1;
            answer= answer2;
            choice= 0;
        elseif choice == 2;
            choice = 0;
        else
            answer= option;
        endif
    else
        if choice == 2;
            choice= 0;
        endif
    clc, %clean shell
    disp("________________________________________________________"), 
    disp("Radiofrequency Circuits and Systems - Useful expressions"), 
    disp("________________________________________________________"), 
    disp(" "), 
    disp("01.-  Input impedance and propagation constant of an "), 
    disp("      infinite-length transmission line"), 
    disp("02.-  T-parameter matrix of a transmission line"), 
    disp("03.-  Lossless transmission line"), 
    disp("04.-  Lossless transmission line in sinusoidal steady state"), 
    disp("05.-  Reflection coefficient"), 
    disp("06.-  Evolution of Ro along the line"), 
    disp("07.-  Voltage wave incident on the load in the general problem"), 
    disp("08.-  Voltage at a distance l respect to the load"), 
    disp("09.-  Current (towards the load) at a distance l respect"), 
    disp("      to the load"), 
    disp("10.-  Voltage standing wave ratio"), 
    disp("11.-  Lossly transmission line"), 
    disp("12.-  Attenuation"), 
    disp("13.-  Power transferred to the load"), 
    disp("14.-  Power delivered by the generator"), 
    disp("15.-  Rectangular wave guide"), 
    disp("16.-  Power density at coordinates with respect to the"), 
    disp("      transmitter"), 
    disp("17.-  Effective area"), 
    disp("18.-  Friis transmission equation"), 
    disp("19.-  Array of two antennas with phase shift"), 
    disp("20.-  Gain of a lambda/2 dipole (arm length = lambda/4)"), 
    disp(" "),
    disp("________________________________________________________"), 
    disp("Another useful expressions"),  
    disp("________________________________________________________"),
    disp("21.-  Trigonometrical arctangent"), 
    disp("22.-  L.T. length in C.C."), 
    disp("23.-  L.T. length in O.C."), 
    disp("24.-  _ViL_ = V2⁻ * endarreriment"),

    disp("30.-  alpha = Attenuation / 20log(e¹)"), 
    disp("31.-  mW a dBm ---> 10log(mW)"),
    disp("32.-  dBm a mW ---> 10^(dBm/10)\n"),

    answer= input("What's your selection? "); 
    endif;

    switch answer
    case 1, %Input Z and propagation constant of an infinite-length T.L.
        Z= input("Z (or L)? "); 
        Y= input("Y (or C)? "); 
        Z0= sqrt(Z/Y); 
        G= sqrt(Z*Y);
        if Z0 == 0
            printf("gamma= %f\n", Z0, G);
        else
            printf("Z0 = %f\ngamma= %f\n", Z0, G);
        endif
    case 2, %T-parameter matrix of a transmission line
        disp("\n      |   s11= cosh(theta)        s12= Z0*sinh(theta)   |"),
        disp(  "  T = |                                                 |"),
        disp(  "      |   s21= 1/sinh(theta)      s22= cosh(theta)      |\n"),
    case 3, %Lossless transmission line
        disp("\ttheta = gamma*l = Tau*s = (l/Vp)*s,   Vp = 1/sqrt(L*C)"),
        yesno= input("Calculate theta (1) or Vp (2)? "), 
        switch yesno
        case 1,
            l= input("l? ");
            G= input("gamma? ");
            Tau= input("Tau? ");
            s= input("s? ");
            Vp= input("Vp? ");
            if Vp == 0, 
                disp("ERROR: Vp can't be zero\n"),
            else
                O1= G*l;
                printf("theta = gamma * l = %f\n", O1);
                disp(""),
                O2= Tau*s;
                printf("theta = Tau * s = %f\n", O2);
                disp(""),
                O3= (l/Vp)*s;
                printf("theta = (l/Vp)*s = %f\n", O3);
                disp(""),
            endif
        case 2,
            Z= input("Z (or L)? "); 
            Y= input("Y (or C)? "); 
            if sqrt(Z*Y) == 0
                disp("CAUTION: Zero division\n");
            else
                Vp= 1/(sqrt(Z*Y));
                printf("Vp = %f\n", Vp);
            endif
        endswitch            
    case 4, %Lossless transmission line in sinusoidal steady state
        Vp= input("Vp? "); 
        f= input("f? "); 
        if f == 0, 
            disp("ERROR: f can't be zero\n"),
        else
            lambda= Vp/f;
            printf("lambda = %f\n", lambda);
        endif
    case 5, %Reflection coefficient
        disp("\tRo = (Z - Z0)/(Z + Z0),       Z = Z0 * (1 + Ro)/(1 - Ro)"),
        yesno= input("Calculate Ro (1) or Z (2)? "), 
        switch yesno
        case 1,
            Z= input("Z? "); 
            Z0= input("Z0? "); 
            if (Z+Z0) == 0,
                disp("CAUTION: Zero division\n");
            else
                Ro= (Z-Z0)/(Z+Z0);
                printf("Ro = %f\n", Ro);
            endif
        case 2,
            Z0= input("Z0? "); 
            Ro= input("Ro? "); 
            if (1-Ro) == 0,
                disp("CAUTION: Zero division\n");
            else
                Z= (1+Ro)/(1-Ro);
                printf("Z = %f\n", Z);
            endif
        otherwise
            error ("invalid value");
        endswitch
    case 6, %Evolution of Ro along the line
        RoL= input("RoL? ");
        alpha= input("alpha? ");
        l= input("l? ");
        theta= alpha*l;
        Ro= RoL*exp(-2*alpha*l);
        printf("INFO VALUE: theta = alpha * l = %f\n", theta);
        printf("Ro = %f\n", Ro);
    case 7, %Voltage wave incident on the load in the general problem
        RoG= input("RoG? ");
        RoL= input("RoL? ");
        alpha= input("alpha? ");
        l= input("l? ");
        theta= alpha*l;
        if (1-RoG*RoL*exp(-2*alpha*l)) == 0,
            disp("CAUTION: Zero division\n");
        else
            VG= input("VG? ");
            ViL= (1/2)*(1-RoG)*(exp(-alpha*l)/(1-RoG*RoL*exp(-2*alpha*l)))*VG;
            printf("INFO VALUE: theta = alpha * l = %f\n", theta);
            printf("ViL = %f\n", ViL);
        endif
    case 8, %Voltage at a distance l respect to the load
        RoL= input("RoL? ");
        alpha= input("alpha? ");
        l= input("l? ");
        theta= alpha*l;
        ViL= input("ViL? ");
        V= ViL*(exp(theta)+RoL*exp(-theta));
        printf("INFO VALUE: theta = alpha * l = %f\n", theta);
        printf("V = %f\n", V);
    case 9, %Current (towards the load) at a distance l respect to the load
        ViL= input("ViL? ");
        Z0= input("Z0? "); 
        if Z0 == 0, 
            disp("ERROR: Z0 can't be zero\n"),
        else
            RoL= input("RoL? ");
            alpha= input("alpha? ");
            l= input("l? ");
            theta= alpha*l;
            I= (ViL/Z0)*(exp(theta)-RoL*exp(-theta));
            printf("INFO VALUE: theta = alpha * l = %f\n", theta);
            printf("I = %f\n", I);
        endif
    case 10, %Voltage standing wave ratio
        Ro= input("Ro? ");
        if (1-abs(Ro)) == 0,
            disp("CAUTION: Zero division\n");
        else
            ROE= (1+abs(Ro))/(1-abs(Ro));
            printf("VSWR = ROE = %f\n", ROE);
        endif
    case 11, %Lossly transmission line
        disp("\ttheta = gamma*l = alpha*l = Tau*s"),
        G= input("gamma? ");
        l= input("l? ");
        alpha= input("alpha? ");
        O1= alpha*l;
        O2= G*l;
        Tau= input("Tau? ");
        s= input("s? ");
        O3= alpha*l+Tau*s;
        printf("INFO VALUE: alpha * l = %f\n", O1);
        disp(""),
        printf("theta = gamma * l = %f\n", O2);
        disp(""),
        printf("theta = alpha * l + Tau * s = %f\n", O3);
        disp(""),
    case 12, %Attenuation
        alpha= input("alpha? ");
        A= alpha*20*log10(exp(1));
        printf("Attenuation (dB/m) = %f\n", A);
    case 13, %Power transferred to the load
        ViL= input("ViL? ");
        Z0= input("Z0? ");
        if Z0 == 0,
            disp("ERROR: Z0 can't be zero\n"),
        else
            RoL= input("RoL? ");
            PL= (ViL^2/(2*Z0))*(1-RoL^2);
            printf("PL = %f [W]\n", PL);
        endif
    case 14, %Power delivered by the generator 
        ViL= input("ViL? ");
        Z0= input("Z0? ");
        if Z0 == 0,
            disp("ERROR: Z0 can't be zero\n"),
        else
            RoL= input("RoL? ");
            alpha= input("alpha? ");
            l= input("l? ");
            PG= (ViL^2/(2*Z0))*(exp(2*alpha*l)-RoL^2*exp(-2*alpha*l));
            printf("PG = %f\n", PG);
        endif
    case 15, %Rectangular wave guide
        lambda= input("lambda? ");
        alpha= input("alpha? ");
        l= input("l? ");
        theta= alpha*l;
        printf("INFO VALUE: theta = alpha * l = %f\n", theta);
        if cos(theta) == 0
            disp("CAUTION: Zero divison\n");
        else
            a= lambda/(2*cos(theta));
            printf("a = %f\n", a);
        endif
    case 16, %Power density at coordinates with respect to the transmitter
        disp("Revisar aquesta fórmula..."),
        PT= input("PT? ");
        r= input("r? ");
        if r == 0
            disp("ERROR: r can't be zero\n");
        else
            GT= input("GT? ");
            alpha= input("alpha? ");
            l= input("l? ");
            theta= alpha*l;
            printf("INFO VALUE: theta = alpha * l = %f\n", theta);
            P= PT/(4*pi*r^2)*GT;
            printf("P = %f\n", P);
        endif
    case 17, %Effective area
        lambda= input("lambda? ");
        G= input("G? ");
        A= (lambda^2/(4*pi))*G
        printf("A = %f\n", A);
    case 18, %Friis transmission equation
        r= input("r? ");
        if r == 0
            disp("ERROR: r can't be zero\n");
        else
            PT= input("PT? ");
            GT= input("GT? ");
            GR= input("GR? ");
            lambda= input("lambda? ");
            PR= PT*GT*GR*(lambda/(4*pi*r)^2);
            printf("PR = %d (W)\n", PR);
            printf("PR = %d (mW)\n", PR*1000);
        endif
    case 19, %Array of two antennas with phase shift
        alpha= input("alpha? ");
        k= input("k? ");
        d= input("d? ");
        theta= input("theta? ");
        printf("INFO VALUE: theta = alpha * l = %f\n", theta);
        FA= 1+exp(j*alpha)*exp(j*k*d*cos(theta));
        printf("FA = %f\n", FA);
    case 20, %Gain of a lambda/2 dipole (arm length = lambda/4)
        alpha= input("alpha? ");
        l= input("l? ");
        theta= alpha*l;
        printf("INFO VALUE: theta = alpha * l = %f\n", theta);
        G= 1.64*(cos((pi/2)*cos(theta))/sin(theta))^2;
        printf("G = %f\n", G);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%    Another useful expressions      %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 21,
        Vp= input("Vp? "); 
        f= input("f? "); 
        Z0= input("Z0? "); 
        Zin= input("Zin? "); 
        l= (atan(Z0/Zin)*Vp)/(2*pi*f);
        printf("l = %f\n", l);
    case 22, %Curt Circuit
        Z= input("Z (or L)? "); 
        Y= input("Y (or C)? "); 
        f= input("f? "); 
        Zin= input("Zin? "); 
        Vp= 1/sqrt(Z*Y);
        printf("Vp = %f\n", Vp);
        Z0= j*sqrt(Z/Y);
        printf("Z0 = %f\n", Z0);
        Tw= (2*pi*f)/Vp;
        printf("Tw = %f\n", Tw);
        l= (atan(Zin/Z0)/Tw);
        printf("l = %f\n", l);
        l2= l+(Vp/f)/2 %if l < 0
        printf("l2 = %f\n", l2);
    case 23, %Circuit Obert
        Z= input("Z (or L)? "); 
        Y= input("Y (or C)? "); 
        f= input("f? "); 
        Zin= input("Zin? "); 
        Vp= 1/sqrt(Z*Y);
        printf("Vp = %f\n", Vp);
        Z0= -j*sqrt(Z/Y);
        printf("Z0 = %f\n", Z0);
        Tw= (2*pi*f)/Vp;
        printf("Tw = %f\n", Tw);
        l= (atan(Z0/Zin)/Tw);
        printf("l = %f\n", l);
    case 24,
        V2n= input("V2⁻? "); 
        alpha= input("alpha? ");
        l= input("l? ");
        theta= alpha*l;
        ViL_= V2n*exp(-theta);
        printf("_ViL_ = %f\n", ViL_);




    case 30, 
        disp("alpha = Attenuation / 20log(e¹)"), 
        A= input("A? ");
        alpha= A/(20*log10(exp(1)));
        printf("alpha = %f\n", alpha);
    case 31, 
        disp("mW a dBm ---> 10log(mW)"), 
        mW= input("mw? ");
        dBm= 10*log10(mW);
        printf("dBm = %f\n", dBm);
        remain = 0;
    case 32, 
        disp("dBm a mW ---> 10^(dBm/10)"), 
        dBm= input("dBm? ");
        mW= 10^(dBm/10);
        printf("mW = %d\n", mW);


    case 40,
        lambda= input("lambda? ");
        PR= input("PR? "); 
        PT= input("PT? "); 
        GT= input("GT? "); 
        GR= input("GR? "); 
        r= lambda/(4*pi*sqrt(PR/(PT*GT*GR)));
        printf("r = %d\n", r);
        r= sqrt((PT*GT*GR*lambda^2)/((4*pi)^2*PR));
        printf("r = %d\n", r);

    otherwise
        error ("invalid value");
    endswitch

    choice= input("Option:\n\tNumber's menu (1)\tMenu (2)\tExit (3): "), 
    switch choice, 
    case 1, 
        answer2 = input("What's your selection? "); 
        remain = 1;
    case 2,
        choice = 2;
        remain = 1;
    case 3,
        remain = 0;
        disp("Bye!!!\n"),

%{
        yesno= input("Do you want to reset variables (y/n)? ", "s"); %"s" indica string
        switch yesno, 
        case {"Yes" "yes" "YES" "y" "Y" "SI" "Si" "S" "si" "s"}, %"s" indica string
            clear,
            remain = 1,
        case {"No" "no" "NO" "n" "N"}, 
            remain = 1;
        otherwise
            error ("invalid value");
        endswitch, 
%}

    otherwise
      error ("invalid value");
    endswitch

    endwhile
    
endfunction



