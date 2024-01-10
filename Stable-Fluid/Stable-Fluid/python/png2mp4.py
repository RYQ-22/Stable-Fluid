import cv2
import argparse
import os

def images_to_video(image_folder, output_file, fps=30):
    # get images and sort
    images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
    images.sort(key=lambda x: int(x.split('_')[1].split('.')[0]))

    # get shape
    frame = cv2.imread(os.path.join(image_folder, images[0]))
    height, width, layers = frame.shape

    # config video writer
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    video_writer = cv2.VideoWriter(output_file, fourcc, fps, (width, height))

    # write mp4
    for image in images:
        img = cv2.imread(os.path.join(image_folder, image))
        video_writer.write(img)

    # release
    video_writer.release()

if __name__ == "__main__":
    # input
    input_folder = "C:/Users/11862/Desktop/vs_code/Fluid-Simulation/Stable-Fluid/Stable-Fluid/Stable-Fluid/python/png"
    # output
    parser = argparse.ArgumentParser(description='get output_path')
    parser.add_argument('-output_file', type=str, help='Description for output file')
    parser.add_argument('-fps', type=int, help='Description for frame per second', default=30)
    args = parser.parse_args()
    output_file = args.output_file
    fps = args.fps
    
    images_to_video(input_folder, output_file, fps)
