
module Import_Test

using FactCheck, JNeuron, TestSetup

facts() do

    myimport=input(string(dirname(Base.source_path()),"/../examples/data/cell2.asc"));
    @fact length(myimport.sections) --> 198

end


end
