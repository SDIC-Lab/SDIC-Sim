function XCOMdata = readXCOMdata(~,element)
    % DESCRIPTION: 
    % read XCOM data from web database based on element name

    date = "2024/01";
    if strcmp(element,"W")
        date = "2024/09";
    end 
    cross_section = urlread("https://sdiclab.deib.polimi.it/wp-content/uploads/"+date+"/XCOM_"+element+"_1keV_500keV.txt");
    cross_section = splitlines(cross_section);
    cross_section = cross_section(~cellfun('isempty',cross_section)); 
    % split based on tab delimiter 
    cross_section = split(cross_section," ",2);  
    XCOMdata.energy = str2num(cell2mat(cross_section(2:end,1)));
    XCOMdata.compton = str2num(cell2mat(cross_section(2:end,2)));
    XCOMdata.photoelectric = str2num(cell2mat(cross_section(2:end,3)));
    XCOMdata.total = str2num(cell2mat(cross_section(2:end,4)));
end 

