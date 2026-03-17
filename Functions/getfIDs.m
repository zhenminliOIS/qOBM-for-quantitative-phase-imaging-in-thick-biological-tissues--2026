function [fIDs,J] = getfIDs(files,bigEndian)
    fIDs = {};
    f = 0; k = 0; f1=0;
    for n=1:length(files)
        if files(n).bytes == 0; continue; end
        tokens = regexp(files(n).name,'(.*)_(\d+)_(\d+)__X(-?\d+\.?\d*)_Y(-?\d+\.?\d*)__Z(-?\d+\.\d+)_(\d+\.\d+)_(.*).bin','tokens');

        tokens = regexp(files(n).name,'(.*)_(\d+)_(\d+)__X(-?\d+\.?\d*)_Y(-?\d+\.?\d*)__Z(-?\d+\.\d+)_(\d+\.\d+)_(.*).bin','tokens');
        if ~isempty(tokens)
            f = f + 1;
            J(f).ID = strcat(tokens{1}{1},'_',tokens{1}{2});
            %J(f).ID = tokens{1}{1};
            J(f).N = (str2num(tokens{1}{3}));%.*4)+str2num(tokens{1}{3});
            J(f).X = str2num(tokens{1}{4});
            J(f).Y = str2num(tokens{1}{5});
            J(f).Z = str2num(tokens{1}{6});
            J(f).qual = str2num(tokens{1}{7});
            J(f).mlb = ~bigEndian;
            J(f).time = datenum(tokens{1}{8},'yy-mm-dd_HH-MM-SS');
            J(f).fname = files(n).name;
        else
            tokens = regexp(files(n).name,'(.*)_(\d+)__X(-?\d+\.?\d*)_Y(-?\d+\.?\d*)__Z(-?\d+\.\d+)_(\d+\.\d+)_(.*).bin','tokens');
            if ~isempty(tokens)
                f = f + 1;
                J(f).ID = tokens{1}{1};
                J(f).N = str2num(tokens{1}{2});
                J(f).X = str2num(tokens{1}{3});
                J(f).Y = str2num(tokens{1}{4});
                J(f).Z = str2num(tokens{1}{5});
                J(f).qual = str2num(tokens{1}{6});
                J(f).mlb = ~bigEndian;
                J(f).time = datenum(tokens{1}{7},'yy-mm-dd_HH-MM-SS');
                J(f).fname = files(n).name;
            else
                % Check if its a cap first
                tokens = regexp(files(n).name,'(.*)_capNum(\d+)_(\d+)_time(\d+\.?\d*).bin','tokens');
                if ~isempty(tokens)
                    f = f + 1;
                    J(f).ID = [tokens{1}{1},'_capNum',tokens{1}{2}];
                    J(f).N = str2num(tokens{1}{3});
                    J(f).X = 0;
                    J(f).Y = 0;
                    J(f).Z = 0;
                    J(f).mlb = ~bigEndian;
                    J(f).time = str2num(tokens{1}{4});
                    J(f).fname = files(n).name;
                else
                    tokens = regexp(files(n).name,'(.*)capNum(\d+)_(\d+)_time(\d+).bin','tokens');
                    if ~isempty(tokens)
                        f = f + 1;
                        J(f).ID = strcat(tokens{1}{1},'capNum',tokens{1}{2});
                        J(f).N = str2num(tokens{1}{3});
                        J(f).X = 0;
                        J(f).Y = 0;
                        J(f).Z = 0;
                        J(f).mlb = ~bigEndian;
                        J(f).time = str2num(tokens{1}{4});
                        J(f).fname = files(n).name;
                    else
                        tokens = regexp(files(n).name,'(.*)capNum(\d+)_(\d+)_time(\d+\.?\d*).bin','tokens');
                        if ~isempty(tokens)
                            f = f + 1;
                            J(f).ID = strcat(tokens{1}{1},'capNum',tokens{1}{2});
                            J(f).N = str2num(tokens{1}{3});
                            J(f).X = 0;
                            J(f).Y = 0;
                            J(f).Z = 0;
                            J(f).mlb = ~bigEndian;
                            J(f).time = str2num(tokens{1}{4});
                            J(f).fname = files(n).name;
                        else
                            tokens = regexp(files(n).name,'(.*)_(\d+)_time(\d+\.?\d*).bin','tokens');
                            if ~isempty(tokens)
                                f = f + 1;
                                J(f).ID = tokens{1}{1};
                                J(f).N = str2num(tokens{1}{2});
                                J(f).X = 0;
                                J(f).Y = 0;
                                J(f).Z = 0;
                                J(f).mlb = ~bigEndian;
                                J(f).time = str2num(tokens{1}{3});
                                J(f).fname = files(n).name;
                            else
                                tokens = regexp(files(n).name,'(.*)_(\d+)_time0.bin','tokens');
                                if ~isempty(tokens)
                                    f = f + 1;
                                    J(f).ID = tokens{1}{1};
                                    J(f).N = str2num(tokens{1}{2});
                                    J(f).X = 0;
                                    J(f).Y = 0;
                                    J(f).Z = 0;
                                    J(f).mlb = ~bigEndian;
                                    J(f).time = 0;
                                    J(f).fname = files(n).name;
                                else
                                    tokens = regexp(files(n).name,'(.*)_(\d+)__X(-?\d+\.?\d*)_Y(-?\d+\.?\d*)__Z(-?\d+\.\d+)_(.*).bin','tokens');
                                    if ~isempty(tokens)
                                        f = f + 1;
                                        J(f).ID = tokens{1}{1};
                                        J(f).N = str2num(tokens{1}{2});
                                        J(f).X = str2num(tokens{1}{3});
                                        J(f).Y = str2num(tokens{1}{4});
                                        J(f).Z = str2num(tokens{1}{5});
                                        J(f).mlb = ~bigEndian;
                                        J(f).time = datenum(tokens{1}{6},'yy-mm-dd_HH-MM-SS');
                                        J(f).fname = files(n).name;
                                    else
                                        tokens = regexp(files(n).name,'(.*)_(\d+)__X(-?\d+\.?\d*)_Y(-?\d+\.?\d*)_Z(-?\d+\.\d+)_time(\d+\.?\d*).bin','tokens');
                                        if ~isempty(tokens)
                                            f = f + 1;
                                            J(f).ID = tokens{1}{1};
                                            J(f).N = str2num(tokens{1}{2});
                                            J(f).X = str2num(tokens{1}{3});
                                            J(f).Y = str2num(tokens{1}{4});
                                            J(f).Z = str2num(tokens{1}{5});
                                            J(f).time = str2num(tokens{1}{6});
                                            J(f).mlb = ~bigEndian;
                                            J(f).fname = files(n).name;
                                        else
                                            tokens = regexp(files(n).name,'(.*)_(\d+)__X(-?\d+\.?\d*)_Y(-?\d+\.?\d*)_Z(-?\d+\.\d+)_time(\d+\.?\d*).bin','tokens');
                                            if ~isempty(tokens)
                                                f = f + 1;
                                                J(f).ID = tokens{1}{1};
                                                J(f).N = str2num(tokens{1}{2});
                                                J(f).X = str2num(tokens{1}{3});
                                                J(f).Y = str2num(tokens{1}{4});
                                                J(f).Z = str2num(tokens{1}{5});
                                                J(f).time = str2num(tokens{1}{6});
                                                J(f).mlb = ~bigEndian;
                                                J(f).fname = files(n).name;
                                            else
                                                % For not fully labeled data from early matlab aq app stages
                                                tokens = regexp(files(n).name,'(.*)_(\d+).bin','tokens');
                                                if ~isempty(tokens)
                                                    f = f + 1;
                                                    J(f).ID = tokens{1}{1};
                                                    J(f).N = str2num(tokens{1}{2});
                                                    J(f).X = 0;
                                                    J(f).Y = 0;
                                                    J(f).Z = 0;
                                                    J(f).mlb = ~bigEndian;
                                                    J(f).time = 0;
                                                    J(f).fname = files(n).name;
                                                else
                                                    % For old data
                                                    tokens = regexp(files(n).name,'(.*)_(\d+)__X(-?\d+\.?\d*)_Y(-?\d+\.?\d*)_(.*).bin','tokens');
                                                    if ~isempty(tokens)
                                                        f = f + 1;
                                                        J(f).ID = tokens{1}{1};
                                                        J(f).N = str2num(tokens{1}{2});
                                                        J(f).X = 0;
                                                        J(f).Y = 0;
                                                        J(f).Z = 0;
                                                        J(f).mlb = ~bigEndian;
                                                        J(f).time = 0;
                                                        J(f).fname = files(n).name;
                                                    else


                                                        tokens = regexp(files(n).name,'(.*)_angles.txt','tokens');
                                                        if ~isempty(tokens)
                                                            f1 = f1 + 1;
                                                            fIDs{f1} = tokens{1}{1};
                                                        else
                                                            tokens = regexp(files(n).name,'(.*)_(\d+)__X(-?\d+\.?\d*)_Y(-?\d+\.?\d*)_Z(-?\d+\.\d+)_(.*).bin','tokens');
                                                            if ~isempty(tokens)
                                                                f = f + 1;
                                                                J(f).ID = tokens{1}{1};
                                                                J(f).N = str2num(tokens{1}{2});
                                                                J(f).X = str2num(tokens{1}{3});
                                                                J(f).Y = str2num(tokens{1}{4});
                                                                J(f).Z = str2num(tokens{1}{5});
                                                                %                                             J(f).qual = str2num(tokens{1}{6});
                                                                J(f).mlb = ~bigEndian;
                                                                J(f).time = datenum(tokens{1}{6},'yy-mm-dd_HH-MM-SS');
                                                                J(f).fname = files(n).name;
                                                            else %% probe 2+ 
                                                                tokens = regexp(files(n).name,'(.*)_(\d+)__ETLmA(.*)!_(.*).bin','tokens');
                                                                if ~isempty(tokens)
                                                                    f = f + 1;
                                                                    J(f).ID = tokens{1}{1};
                                                                    J(f).N = str2num(tokens{1}{2});
                                                                    J(f).Z = str2num(tokens{1}{3});
                                                                    J(f).X = 0; % dummy
                                                                    J(f).Y = 0; % dummy
                                                                    %                                             J(f).qual = str2num(tokens{1}{6});
                                                                    J(f).mlb = ~bigEndian;
                                                                    J(f).time = datenum(tokens{1}{4},'yy-mm-dd_HH-MM-SS');
                                                                    J(f).fname = files(n).name;
                                                                    etlTag = true;
                                                                else
                                                                    disp([files(n).name,' cannot be tokenize correctly\n'])
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end