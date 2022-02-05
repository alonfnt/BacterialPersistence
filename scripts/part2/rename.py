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
        _, dst = src.split('-first')
        dst = f"{dst}-min"
        dst = f"{folder}/mutation/{dst}"
         
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