# Pythono3 code to rename multiple
# files in a directory or folder
 
# importing os module
import os
 
# Function to rename multiple files
def main():
   
    folder = "../../data/model3/"
    for filename in os.listdir(folder+"rename/"):
        src = f"{filename}"  # foldername/filename, if .py file is outside folder
        print(src)

        dst = src.split('0.')
        dst = dst[0] + dst[1]
        # dst = src.split('-T0')
        # dst = dst[0] + "0-T0" + dst[1]

        # if len(src.split('l-')) == 2:
        #     _, dst = src.split('l-')
        #     dst = f"{dst}-optimal"

        # if len(src.split('t-')) == 2:
        #     _, dst = src.split('t-')
        #     dst = f"{dst}-min"
        
        # if len(src.split('x-')) == 2:
        #     _, dst = src.split('x-')
        #     dst = f"{dst}-max"

        dst = f"{folder}/{dst}"

         
        # rename() function will
        # rename all the files
        os.rename(folder+"rename/"+src, dst)

    for filename in os.listdir(folder+"rename/"):
        src = f"{filename}"  # foldername/filename, if .py file is outside folder
        os.remove(folder+"rename/"+src)
 
# Driver Code
if __name__ == '__main__':
     
    # Calling main() function
    main()