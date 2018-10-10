    function [] = plot_dubins(a, b, q1, q2, x, y, H, gamma1, gamma3, sgn1, sgn3, type, s1, s3)
    
        global LSL LSR RSL RSR RLR LRL TSL TSR LST RST TST LSEG SSEG RSEG TSEG INF TYPENAME DIRNAME DIRDATA;

        dubins_const;
        
        i = type + 1;
        fprintf(1, 'type=%d i=%d H=%f gamma1=%f gamma3=%f sgn1=%d sgn3=%d s1=%f s3=%f\n', type, i, H, gamma1, gamma3, sgn1, sgn3, s1, s3);
        
        px=[a(1), q1(1), q2(1), b(1)];
        py=[a(2), q1(2), q2(2), b(2)];
        
        xmin = min(px); 
        ymin = min(py); 
        xmax = max(px); 
        ymax = max(py); 
    
        xmin1 = min(x); 
        ymin1 = min(y); 
        xmax1 = max(x); 
        ymax1 = max(y); 

        if (xmin1 < xmin) xmin = xmin1; end;
        if (ymin1 < ymin) ymin = ymin1; end;
        if (xmax1 > xmax) xmax = xmax1; end;
        if (ymax1 > ymax) ymax = ymax1; end;
        
        dx = b(1) - a(1);
        dy = b(2) - a(2);
        D = sqrt(dx*dx + dy*dy);
        theta = wrapTo2Pi(atan2(dy, dx));
        
        %xmin = -150.0;
        %xmax = 300.0;
        %ymin = -110.0;
        %ymax = 120.0;
    
        H2 = D*0.16;
        xmin = xmin - H2;
        xmax = xmax + H2;
        ymin = ymin - H2;
        ymax = ymax + H2;
        fprintf(1, 'xmin=%f xmax=%f ymin=%f ymax=%f\n', xmin, xmax, ymin, ymax);
                    
        plot(px, py, 'rs');
        axis([xmin xmax ymin ymax]);
        axis equal;    
        %set(gca,'XTick',-pi:pi/2:pi)
        %set(gca, 'ylim', [-2 5], 'xlim', [-2 5]);
        hold all;
        
        xLimits = get(gca,'XLim');  %# Get the range of the x axis
        yLimits = get(gca,'YLim');  %# Get the range of the y axis
        fprintf(1, 'xLimits: %f %f\n', xLimits(1), xLimits(2));
        fprintf(1, 'yLimits: %f %f\n', yLimits(1), yLimits(2));
           
        plot(x, y);
    
        k=D*0.08;
        text(xLimits(1)+k, yLimits(2)-k, TYPENAME(i,:), 'FontSize', 12, 'fontWeight','bold');
    
        gamma1_ = gamma1 + theta;
        gamma3_ = gamma3 + theta;
        of = 4*H;
           
        x2(1) = a(1);
        x2(2) = a(2);
        y2(1) = a(1) + cos(a(3))*k;
        y2(2) = a(2) + sin(a(3))*k;
        arrow(x2, y2, 'Length', 8, 'BaseAngle', 10, 'TipAngle', 20);
        %arrow FIXLIMITS;
        %ARROW(Start,Stop,Length,BaseAngle,TipAngle,Width,Page,CrossDir)
        text(y2(1)+of, y2(2)+of,' S', 'FontSize',12, 'fontWeight','bold');
        
        x2(1) = b(1);
        x2(2) = b(2);
        y2(1) = b(1) + cos(b(3))*k;
        y2(2) = b(2) + sin(b(3))*k;
        arrow(x2, y2, 'Length', 8, 'BaseAngle', 10, 'TipAngle', 20);
        %arrow FIXLIMITS;
        text(y2(1)+of, y2(2)+of,' G', 'FontSize',12, 'fontWeight','bold');
    
        xlabel('x', 'FontSize',12, 'fontWeight','bold');
        ylabel('y', 'FontSize',12, 'fontWeight','bold');
    
        if (DIRDATA(i, 1) == TSEG)    
            plot_turn_in_place(k, H, sgn1, gamma1_, a(1), a(2), a(3), s1)
        end;
    
        if (DIRDATA(i, 3) == TSEG)
            plot_turn_in_place(k, H, sgn3, gamma3_, b(1), b(2), gamma3_, s3)
        end;
       
        hold off;
        grid on;
    end