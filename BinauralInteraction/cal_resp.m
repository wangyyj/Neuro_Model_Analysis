function cal_resp(hObject,eventdata,handles)
global g_speaker
g_speaker.cal_response=[];
SR_AD=g_speaker.SR_AD2; 

if hObject==handles.cal_resp_start
    g_speaker.play=1; 
    set(handles.cal_resp_start,'visible','off'); 
%     set(handles.pb_stop,'visible','on'); 

    if g_speaker.type==1  % golay code
        
        %%%%%%%%%%%% Convolve stim with inverse filter
            y0 = y;
            clear y;
            y = 
            y = conv(y,hz_inv_mp);
        
        
        
        % a and b sequence responses
     	if g_speaker.demo == 0
            for stim_num = 1:2
              	start_tdt(stim_num);
                response(stim_num,:) = speaker_engine*10/32767;	
                                                    % raw AD data in int16, converted into volt
                                                    % 32767 = 2^15-1, +-10V danamic range on RX6
            end;
        else
                response = speaker_simulated_engine;
        end;
        g_speaker.response = response;   	
       
        % Golay Code Algorithm; only use 2^N=L points for FFT 
        L = g_speaker.Lsequence; 
        L_temp = round(log2(L)); 
        [a,b] = generate_golay(L_temp);             % pulse voltage is unit 1

        % remove offset samples and voltage
        response_a = g_speaker.response(1,g_speaker.sounddelay+1:g_speaker.sounddelay+L);           
        response_b = g_speaker.response(2,g_speaker.sounddelay+1:g_speaker.sounddelay+L);
        response_a = response_a-g_speaker.response(1,1); 
        response_b = response_b-g_speaker.response(2,1); 
        response =  [response_a; response_b];

        % get Golay FFT response in frequency domain, this line is taken from EY's code
        system_response = (fft(response_a).*conj(fft(a))+fft(response_b).*conj(fft(b)))/(2*L);

        % test reverse filter (revere filter and minimum-phase reverse filter)
        [g_speaker.h_direct, g_speaker.hz_inv,g_speaker.hz_inv_mp] = ...
            h_reverse(system_response, SR_AD);

        % prepare a compatible amplitude system_response for following shared analysis
        system_response =  fft(response')';       

    else % click, tone, BP-noise, user
        if g_speaker.demo == 0
            start_tdt(1);
            %%%%%%%%%%%% Convolve stim with inverse filter

            
            
            
            response = speaker_engine*10/32767;     % raw AD data in int16, converted into volt
                                                    % 32767 = 2^15-1, +-10V danamic range on RX6
        else
            response = speaker_simulated_engine;
        end;
       	g_speaker.response = response;

        % calculate FFT using AC component (not in use)
        if g_speaker.type==2
            L = g_speaker.Lsequence + round(g_speaker.SR_AD2*g_speaker.post_stim);
        elseif g_speaker.type==3 || g_speaker.type==4
            L = round(g_speaker.SR_AD2*g_speaker.tone_dur);
        else
            L = round(g_speaker.SR_AD2*g_speaker.user_dur);
        end;

        response = g_speaker.response(:,g_speaker.sounddelay+1:g_speaker.sounddelay+L)...
            -g_speaker.response(1);
        system_response = fft(response')';
    end % end if gloay
      
    % Power of the whole response, temporal / spectral
    p_response = response*1000;                    	% mV, convert V to mV
    p_response = p_response/g_speaker.filter_gain; 	% mV, signal before pre-amp filter
    p_response = p_response/g_speaker.mVPa;        	% Pascal, acoutic amplitude
    p_response = sqrt(mean(p_response.^2'));     	% Pascal rms from temporal analysis  
    p_response = sqrt(mean(p_response.^2));         % average across    a&b for golay
                                                    %                   repetitions for others
                                                    
    P_response = system_response/L;                 % V/freq bin, nomalized by L 
    P_response = abs(P_response(:,1:floor(L/2)+1));
	P_response(:,2:ceil(L/2)) = P_response(:,2:ceil(L/2))*2;             	
                                                    % V/freq bin, combined mirrored amplitude
    P_response = P_response*1000;                   % mV/freq bin, convert V to mV                  
    P_response = P_response/g_speaker.filter_gain;	% mV/freq bin, signal before pre-amp filter
    P_response = P_response/g_speaker.mVPa;      	% Pascal/freq bin, on microphone
    P_response(:,2:ceil(L/2)) = P_response(:,2:ceil(L/2))/sqrt(2);
                                                    % Pascal rms/freq bin, on microphone
    P_response = sqrt(mean(P_response.^2));         
                                                    
    power_response = 20*log10(P_response/1) + g_speaker.refSPL;
                                                    % dB SPL, abosolute sound level, (1Pascal rms = 94dB SPL)
   	power_response = max(power_response, -40+zeros(size(power_response,1), size(power_response,2)));    
                                                    % get rid of -Inf levels, set baseline to be -40dB SPL 
    g_speaker.power_response = power_response;      % average across    a&b for golay
                                                    %                   repetitions for others

	% display numbers
    disp(['Total Pascal(rms) from temporal analysis is ', num2str(p_response),','])   
    disp(['which equals to ', num2str(20*log10(p_response)+g_speaker.refSPL),' dB SPL'])
    disp(['Total Pascal(rms) from spectral analysis is ', num2str(sqrt(sum(P_response.^2,2)')),','])   
    disp(['which equals to ', num2str(10*log10(sum(P_response.^2,2)')+g_speaker.refSPL),' dB SPL'])
    
    if  g_speaker.type==1 %|| g_speaker.type==2 %% click or golay 
        % Power of the direct response of golay h(101:356) points 
        % change the length in the denominator 01/09/2009 YZ; It is wrong to use L to calculate H_direct
        H=fft(g_speaker.h_direct); 
        L_direct=length(g_speaker.h_direct); 
        power_response2=2*abs(H).^2/L_direct/SR_AD; %% only 1:L_direct/2 is useful
        power_response2(1)=power_response2(1)/2; 
        % adjust for the mic sensitivity and analog filter gain
        power_response2=10*log10(power_response2)+ g_speaker.refSPL - g_speaker.dBV ;

        % adjust for analog filter gain and DC block gain
        power_response2=power_response2 - 20*log10(g_speaker.filter_gain);
        g_speaker.power_response_direct=power_response2; 
    end 
      
   % plot response data
   y_in = g_speaker.y(1,:); 
   y_out = g_speaker.response(end,:); %/max(abs(g_speaker.ad2data(end,:)));  %/(2^15); %*10^(g_speaker.att/20); 
   showit(y_out,power_response,handles,SR_AD); 

   % Show Power, and Freq bin, but not THD yet
   [g_speaker.peak_power, f_index]=max(power_response(10:end)); % avoid DC component 
   f_index = f_index+9;
   power_total_str = sprintf('%6.2f',20*log10(p_response) + g_speaker.refSPL);
   power_peak_str =  sprintf('%6.2f',(g_speaker.peak_power));
   power_5pts_str =  sprintf('%6.2f',10*log10(sum(P_response(max(f_index-2,1):min(f_index+2,L)).^2))+g_speaker.refSPL);  
   
   if g_speaker.type ==3 % tone
      f_max=f_index*SR_AD/L;  
      if abs(f_max-g_speaker.CF)<=g_speaker.CF*0.02; % within 2% CF
        set(handles.power_edit,'string',...
           [power_peak_str,' / ',power_5pts_str,' / ',power_total_str,' (dB SPL)']); 
      else
        set(handles.power_edit,'string','NA');
      end
   else  % click, golay, BP-noise
        set(handles.power_edit,'string',...
           [power_peak_str,'/',power_5pts_str,'/',power_total_str,' (dB SPL)']);
   end     
        set(handles.FreqBin_edit,'string',...
           strcat(num2str(SR_AD/L,'%.2f'), ' (Hz)'));   
   
   %% update gui 
   set(handles.cal_resp_start,'visible','on'); 
%    set(handles.pb_stop,'visible','off'); 
   infor=strcat(g_speaker.type_string(g_speaker.type),' is done!'); 
   set(handles.code_dir,'string',infor); 
   
elseif hObject==handles.pb_stop
    g_speaker.play=0;
    set(handles.pb_start,'visible','on'); 
    set(handles.pb_stop,'visible','off'); 
    if g_speaker.demo == 0
      stop_tdt;
    end     
    infor=strcat(g_speaker.type_string(g_speaker.type),' is not completed!'); 
    set(handles.code_dir,'string',infor); 
end