function paramS = MibiReadRunXml(fileName,pointNumber)
% function params = MibiReadRunXml(fileName,PointNumber)
% function reads run parameters from xml for point number PointNumber

paramNames= {'TimeResolution', 'MassGain', 'MassOffset','XSize', 'YSize', 'PointName'};
params= cell(size(paramNames));

% read xml to string
textXML = fileread(fileName);
% find parameters
for i=1:length(paramNames)
    pattern=[paramNames{i},'="([\+-\w.]+)"\>'];
    [matchExp,tok,ext]= regexp(textXML, pattern, 'match','tokens','tokenExtents');
    [params{i}, status]=str2num(tok{pointNumber}{1});
    if status==0
        params{i}=tok{pointNumber}{1};
    end
end
paramS = cell2struct(params', paramNames');

% % get the xpath mechanism into the workspace
% import javax.xml.xpath.*
% factory = XPathFactory.newInstance;
% xpath = factory.newXPath;
% 
% % read file
% docNode = xmlread(fileName);
% 
% % compile and evaluate the XPath Expression
% expression = xpath.compile('//Depth_Profile');
% DepthProfileNode = expression.evaluate(docNode, XPathConstants.NODE);
% DepthProfileAtt = DepthProfileNode.getAttributes;
% 

% docRoot = xDoc.getDocumentElement;
% 
% phoneNumber = friendlyInfo.getElementsByTagName('PhoneNumber').item(0).getTextContent
% 
% 
% s=xml2struct(fileName);
% allDepthProfiles = docNode.getElementsByTagName('Depth_Profile');
% for k = 0:allDepthProfiles.getLength-1
%    thisListitem = allDepthProfiles.item(k);
%    
%    % Get the label element. In this file, each
%    % listitem contains only one label.
%    thisList = thisListitem.getElementsByTagName('TimeResolution');
%    thisElement = thisList.item(0);
% 
%    % Check whether this is the label you want.
%    % The text is in the first child node.
%    if strcmp(thisElement.getFirstChild.getData, findLabel)
%        thisList = thisListitem.getElementsByTagName('callback');
%        thisElement = thisList.item(0);
%        findCbk = char(thisElement.getFirstChild.getData);
%        break;
%    end
%    
% end
% Display the final results:
% if ~isempty(findCbk)
%     msg = sprintf('Item "%s" has a callback of "%s."',...
%                   findLabel, findCbk);
% else
%     msg = sprintf('Did not find the "%s" item.', findLabel);
% end
% disp(msg);