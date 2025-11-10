function out = compare_files(actualFile, expectedFile, dontCleanUp)
    arguments
        actualFile (1, 1) string
        expectedFile (1, 1) string
        dontCleanUp (1, 1) double = 0
    end
    actual = fileread(actualFile);
    expected = fileread(expectedFile);
    out = strcmp(actual, expected);
    if (dontCleanUp == 0)
        delete(actualFile);
    end
end
