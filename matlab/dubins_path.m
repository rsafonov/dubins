
function [] = dubins_path()

    global LSL LSR RSL RSR RLR LRL LSEG SSEG RSEG INF;
    
    LSL = 1;
    LSR = 2;
    RSL = 3;
    RSR = 4;
    RLR = 5;
    LRL = 6;
    
    LSEG = 0;
    SSEG = 1;
    RSEG = 2;
    
    INF=100000.0;
    
    %DIRDATA = [ 0, 1, 0; 0, 1, 2; 2, 1, 0; 2, 1, 2; 2, 0, 2; 0, 2, 0;];
    DIRDATA = [LSEG, SSEG, LSEG; LSEG, SSEG, RSEG; RSEG, SSEG, LSEG; RSEG, SSEG, RSEG; RSEG, LSEG, RSEG; LSEG, RSEG, LSEG;];
       
    %a = [6.0, 2.0, 2.722];
    %b = [14.0, 1.0, 0.5];
    a = [0.0, 0.0, 0.0];
    b = [6.0, 1.0, 3.142];
    
    %nn=n*rho;
    %h = M_2PI/(nn);

    rho = 1.0;
    %H=0.1;
    H=0.052360;
    fprintf(1, 'a: x=%f y=%f theta=%f\n', a(1), a(2), a(3)); 
    fprintf(1, 'b: x=%f y=%f theta=%f\n', b(1), b(2), b(3)); 
    fprintf(1, 'rho=%f H=%f\n', rho, H); 
        
    t = [INF, INF, INF, INF, INF, INF];
    p = [INF, INF, INF, INF, INF, INF];
    q = [INF, INF, INF, INF, INF, INF];
    cost = [INF, INF, INF, INF, INF, INF];
        
    dx = b(1) - a(1);
    dy = b(2) - a(2);
    D = sqrt(dx * dx + dy * dy);
    d = D/rho;
    fprintf(1, 'dx=%f dy=%f D=%f d=%f\n', dx, dy, D, d);    
    
    theta = wrapTo2Pi(atan2(dy, dx));
    alpha = wrapTo2Pi(a(3) - theta);
    beta  = wrapTo2Pi(b(3) - theta);
    fprintf(1, 'theta=%f alpha=%f beta=%f\n', theta, alpha, beta);

    sa = sin(alpha);            
    sb = sin(beta);             
    ca = cos(alpha);           
    cb = cos(beta);             
    cab = cos(alpha - beta);

    %LSL
    tmp0 = d+sa-sb;
    tmp2 = 2 + (d*d) - (2*cab) + (2*d*(sa - sb));
    if( tmp2 >= 0 && tmp0 >= 0) 
        tmp1 = atan((cb-ca)/tmp0);
        t(LSL) = wrapTo2Pi(-alpha + tmp1);
        p(LSL) = sqrt(tmp2);
        q(LSL) = wrapTo2Pi(beta - tmp1);
        cost(LSL) = t(LSL) + p(LSL) + q(LSL);
    end
    fprintf(1, 'LSL: t=%f p=%f q=%f cost=%f\n', t(LSL), p(LSL), q(LSL), cost(LSL));
        
    %RSR
    tmp0 = d-sa+sb;
    tmp2 = 2 + (d*d) -(2*cab) + (2*d*(sb-sa));
    if( tmp2 >= 0 && tmp0 >= 0) 
        tmp1 = atan( (ca-cb)/tmp0 );
        t(RSR) = wrapTo2Pi(alpha - tmp1);
        p(RSR) = sqrt(tmp2);
        q(RSR) = wrapTo2Pi(-beta + tmp1);
        cost(RSR) = t(RSR) + p(RSR) + q(RSR);        
    end
    fprintf(1, 'RSR: t=%f p=%f q=%f cost=%f\n', t(RSR), p(RSR), q(RSR), cost(RSR));
    
    %LSR
    tmp1 = -2 + (d*d) + (2*cab) + (2*d*(sa+sb));
    if( tmp1 >= 0 ) 
        p(LSR)  = sqrt(tmp1);
        tmp2 = atan((-ca-cb)/(d+sa+sb)) - atan(-2.0/p(LSR));
        t(LSR) = wrapTo2Pi(-alpha + tmp2);
        q(LSR) = wrapTo2Pi(-wrapTo2Pi(beta) + tmp2);
        cost(LSR) = t(LSR) + p(LSR) + q(LSR);
    end    
    fprintf(1, 'LSR: t=%f p=%f q=%f cost=%f\n', t(LSR), p(LSR), q(LSR), cost(LSR));
    
    %RSL
    tmp1 = (d*d) -2 + (2*cab) - (2*d*(sa+sb));
    if( tmp1 > 0 ) 
        p(RSL) = sqrt(tmp1);
        tmp2 = atan((ca+cb)/(d-sa-sb)) - atan(2.0/p(RSL));
        t(RSL) = wrapTo2Pi(alpha - tmp2);
        q(RSL) = wrapTo2Pi(beta - tmp2);
        cost(RSL) = t(RSL) + p(RSL) + q(RSL);
    end
    fprintf(1, 'RSL: t=%f p=%f q=%f cost=%f\n', t(RSL), p(RSL), q(RSL), cost(RSL));
    
    %RLR
    tmp1 = (6.0 - d*d + 2*cab + 2*d*(sa-sb))/8.0;
    if( abs(tmp1) < 1) 
        p(RLR) = acos(tmp1);
        t(RLR) = wrapTo2Pi(alpha - atan2(ca-cb, d-sa+sb) + wrapTo2Pi(p(RLR)/2.0));
        q(RLR) = wrapTo2Pi(alpha - beta - t(RLR) + wrapTo2Pi(p(RLR)));
        cost(RLR) = t(RLR) + p(RLR) + q(RLR);
    end
    fprintf(1, 'RSR: t=%f p=%f q=%f cost=%f\n', t(RLR), p(RLR), q(RLR), cost(RLR));
    
    %LRL
    tmp1 = (6.0 - d*d + 2*cab + 2*d*(-sa + sb))/8.0;
    if( abs(tmp1) < 1) 
        p(LRL) = wrapTo2Pi(acos(tmp1));
        t(LRL) = wrapTo2Pi(-alpha - atan2(ca-cb, d+sa-sb) + p(LRL)/2.0);
        q(LRL) = wrapTo2Pi(wrapTo2Pi(beta) - alpha -t(LRL) + wrapTo2Pi(p(LRL)));
        cost(LRL) = t(LRL) + p(LRL) + q(LRL);
    end
    fprintf(1, 'LRL: t=%f p=%f q=%f cost=%f\n', t(LRL), p(LRL), q(LRL), cost(LRL));

    % Extract the best cost path
    [costMin, iMin] = min(cost); 
    %iMin=4;
    %costMin=cost(iMin);
    length = costMin*rho;
    n=floor(length/H) + 1;
    fprintf(1, 'iMin=%d costMin=%f length=%f n=%d\n', iMin, costMin, length, n);
    fprintf(1, '%d %d %d\n', DIRDATA(iMin, 1), DIRDATA(iMin, 2), DIRDATA(iMin, 3));
        
    a0 = [0.0, 0.0, a(3)];
    fprintf(1, 'a0: %f %f %f\n', a0(1), a0(2), a0(3));
    
    p1 = t(iMin);
    p2 = p(iMin);
    fprintf(1, 'p1=%f p2=%f\n', p1, p2);
    
    q10 = dubins_segment(p1, a0, DIRDATA(iMin, 1));
    q1(1) = q10(1)*rho + a(1);
    q1(2) = q10(2)*rho + a(2);
    q1(3) = wrapTo2Pi(q10(3));
    fprintf(1, 'q1: %f %f %f\n', q1(1), q1(2), q1(3));   
    
    q20 = dubins_segment(p2, q10, DIRDATA(iMin, 2));
    fprintf(1, 'q20: %f %f %f\n', q20(1), q20(2), q20(3));
    q2(1) = q20(1)*rho + a(1);
    q2(2) = q20(2)*rho + a(2);
    q2(3) = wrapTo2Pi(q20(3));
    fprintf(1, 'q2: %f %f %f\n', q2(1), q2(2), q2(3)); 
           
    %x1 = zeros(n);
    %y1 = zeros(n);
    
    tt= 0.0;
    i=1;
    while (tt < length)
        
       tt0 = tt/rho;
       
       if (tt0 < p1) 
           type = DIRDATA(iMin, 1);
           qt = dubins_segment(tt0, a0, type);    
       elseif (tt0 < (p1+p2)) 
           type = DIRDATA(iMin, 2);
           qt = dubins_segment(tt0-p1, q10, type);
       else 
           type = DIRDATA(iMin, 3);
           qt = dubins_segment(tt0-p1-p2, q20, type); 
       end;       
       
       qt(1) = qt(1)*rho + a(1);
       qt(2) = qt(2)*rho + a(2);
       qt(3) = wrapTo2Pi(qt(3));
       
       %fprintf(1, 'i=%d tt=%f x=%f y=%f theta=%f type=%d\n', i, tt, qt(1), qt(2), qt(3), type); 
       x1(i) = qt(1);
       y1(i) = qt(2);
       i = i+1;
       tt  = tt + H;
    end;
    fprintf(1, 'i=%d tt=%f\n', i, tt);

    x=[a(1), q1(1), q2(1), b(1)];
    y=[a(2), q1(2), q2(2), b(2)];
    
    xmin = min(x); 
    ymin = min(y); 
    xmax = max(x); 
    ymax = max(y); 
    
    xmin1 = min(x1); 
    ymin1 = min(y1); 
    xmax1 = max(x1); 
    ymax1 = max(y1); 

    if (xmin1 < xmin) xmin = xmin1; end;
    if (ymin1 < ymin) ymin = ymin1; end;
    if (xmax1 > xmax) xmax = xmax1; end;
    if (ymax1 > ymax) ymax = ymax1; end;
    
    xmin = xmin - H;
    xmax = xmax + H;
    ymin = ymin - H;
    ymax = ymax + H;
    
    fprintf(1, 'xmin=%f xmax=%f ymin=%f ymax=%f\n', xmin, xmax, ymin, ymax);
            
    %axis([-2.0 6.0 -2.0 6.0]);
    plot(x, y, 'rs');
    %axis([xmin xmax ymin ymax]);
    axis equal;    
    %set(gca,'XTick',-pi:pi/2:pi)
    %set(gca, 'ylim', [-2 5], 'xlim', [-2 5]);
    hold all;
    
    plot(x1,y1);
    hold off;
    grid on;
        
    function [qt] = dubins_segment(t, qi, type) 
       if(type == LSEG)
           qt(1) = qi(1) + sin(qi(3)+t) - sin(qi(3));
           qt(2) = qi(2) - cos(qi(3)+t) + cos(qi(3));
           qt(3) = mod2pi(qi(3) + t);
       elseif(type == RSEG)
           qt(1) = qi(1) - sin(qi(3)-t) + sin(qi(3));
           qt(2) = qi(2) + cos(qi(3)-t) - cos(qi(3));
           qt(3) = wrapTo2Pi(qi(3) - t);
       else %type == SSEG 
           qt(1) = qi(1) + cos(qi(3)) * t;
           qt(2) = qi(2) + sin(qi(3)) * t;
           qt(3) = wrapTo2Pi(qi(3));
       end;
    end
end

function [f] = mod2pi(theta)
    f = theta - 2*pi*floor(theta/(2*pi));    
end
