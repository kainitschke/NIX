% start of the NIX toolbox
% no parameters needed
clc;
fprintf('    .----------------. .----------------. .----------------.   \n');
fprintf('   | .--------------. | .--------------. | .--------------. |  \n');
fprintf('   | | ____  _____  | | |     _____    | | |  ____  ____  | |  \n');
fprintf('   | ||_   \\|_   _| | | |    |_   _|   | | | |_  _||_  _| | |  \n');
fprintf('   | |  |   \\ | |   | | |      | |     | | |   \\ \\  / /   | |  \n');
fprintf('   | |  | |\\ \\| |   | | |      | |     | | |    > \\/ <    | |  \n');
fprintf('   | | _| |_\\   |_  | | |     _| |_    | | |  _/  /\\  \\_  | |  \n');
fprintf('   | ||_____|\\____| | | |    |_____|   | | | |____||____| | |  \n');
fprintf('   | |              | | |              | | |              | |  \n');
fprintf('   | ''--------------'' | ''--------------'' | ''--------------'' |  \n');
fprintf('    ''----------------'' ''----------------'' ''----------------''   \n');

%% clear workspace
clear global LUE
k = findall(0,'type','figure');
for i = 1:length(k),
    if isequal(get(k(i),'Name'),'nparLD') | isequal(get(k(i),'Name'),'NIX'),
        delete(k(i));
    end;
end;

%% check for SPM
if isempty(which('spm')),
    uiwait(warndlg(sprintf('SPM was not found. Some subroutines of SPM are needed in the NIX-toolbox. SPM is freeware by the Wellcome Trust Centre for Neuroimaging (www.fil.ion.ucl.ac.uk/spm/).\nYou will now be asked to locate SPM on your system. Without SPM the NIX-toolbox will not work.\nTo avoid this message add SPM permanently to your start paths of Matlab.'),'SPM missing','modal'));
    spm_path = uigetdir(pwd,'Select SPM directory');
    addpath(genpath(spm_path));
end;

%% Start GUI
nix_gui;